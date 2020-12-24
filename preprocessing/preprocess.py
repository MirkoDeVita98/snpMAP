import os
import ntpath
import shutil
from collections import defaultdict

import pandas as pd
from tqdm import tqdm
# import vcf
from time import time 

from utils import file_info, parser, downloader
from constants import interesting_effects, interesting_roles, \
    allowed_alleles, descr_to_label_binary, descr_to_label_multiclass
from io_utils.reading import VCFReader


config = parser.load_config()
cache_dir = config['output']['cache_dir']
dependencies_dir = config['input']['dependencies_dir']


def _hgvs_id(chrom, position, ref, alt):
    return f'chr{chrom}:g.{position}{ref}>{alt}'


def _get_filename(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def get_hgvs_ids(fpath):
    hgvs_fpath = results_fpath(fpath, 'hgvs')
    if not os.path.exists(hgvs_fpath):
        extract_hgvs_ids(fpath)

    return pd.read_csv(hgvs_fpath)['ID'].to_list()


def results_fpath(input_fpath, name, ending='csv', folder=cache_dir):
    input_fname = _get_filename(input_fpath)
    parts = input_fname.split('.')  # MUST satisfy
    base_name = parts[0]
    results_fname = f'{base_name}_{name}{".".join(parts[1:-1])}.{ending}'
    return os.path.join(folder, results_fname)


def shrink_vcf(fpath, mode):
    if mode == 'train':
        print('The train file does not need to be shrunken.')
        return fpath

    shrunken_fpath = results_fpath(fpath, 'shrunken')

    if os.path.exists(shrunken_fpath):
        print('The VCF file has already been shrunken.')
        return shrunken_fpath

    print('Shrinking the VCF file (removing the unnecessary information).')

    flen = file_info.file_length(fpath)
    with open(fpath, 'r') as fr, open(shrunken_fpath, 'w') as fw:
        for line in tqdm(fr, total=flen):
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=AC') \
                or not line.startswith('##INFO=<ID'):
                    fw.write(line)
            elif line.startswith('#CHROM'):
                fw.write(
                    '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n')
            else:
                fw.write(line[:line.find(';')] + '\tGT\n')
                # fw.write(line[:line.find(';') + 1] + \
                #     line[line.find('ANN='): line.find('GT') + 2] + '\n')
    
    print('Shrinking done.')

    return shrunken_fpath


def annotate_vcf(fpath, mode):
    if mode == 'train':
        print('The train file does not need to be annotated.')
        return fpath

    annotated_path = results_fpath(fpath, 'annotated')
    if os.path.exists(annotated_path):
        print('The VCF file has already been annotated with SNP type info.')
        return annotated_path

    if not os.path.exists('dependencies/snpEff/snpEff.jar'):
        downloader.download_SNPEff()

    print('Annotating the VCF file with information about the SNP types.')
    print('>> This may take a few minutes!')

    snpEff_dir = os.path.join(dependencies_dir, 'snpEff')
    
    db_path = os.path.join(snpEff_dir, 'data', 'GRCh37.75')
    if not os.path.exists(db_path):
        print('The database must be extracted first.')
        archive_path = os.path.join(snpEff_dir, 'data', 'gen.tar.gz')
        tmp_path = os.path.join(snpEff_dir, 'data', 'tmp_db')
        shutil.unpack_archive(archive_path, tmp_path)
        shutil.move(os.path.join(tmp_path, 'GRCh37.75'), db_path)
        shutil.rmtree(tmp_path)
        print('Database extracted.')

    start_time = time()
    os.system(
        f'java -Xmx6g -jar dependencies/snpEff/snpEff.jar GRCh37.75 -noHgvs -noLof {fpath} > {annotated_path}')
    end_time = time()
    elapsed_time = end_time - start_time
    print(f'Annotation finished in {int(elapsed_time // 60)} minutes and {int(elapsed_time % 60)} seconds.')

    return annotated_path


def extract_SNPs(annotated_fpath, SNP_type='nsSNP', mode='train'):
    output_fpath = results_fpath(annotated_fpath, SNP_type)

    if os.path.exists(output_fpath):
        print('SNPs have already been extracted.')
        return output_fpath

    print('Extracting SNPs from the annotated file.')

    flen = file_info.file_length(annotated_fpath)
    metalen = file_info.no_info_lines(annotated_fpath)
    no_of_variants = flen - metalen

    snp_count = 0
    effects = defaultdict(int)

    with open(annotated_fpath, 'r') as fr, open(output_fpath, 'w') as fw:
        # for line in tqdm(fr, total=flen):
        #     if line.startswith('#'):
        #         fw.write(line)
        #     elif _SNP_condition(mode, SNP_type, line):
        #         # fw.write(line[:line.find(';') + 1] + '\tGT\n')
        #         # fw.write(line[:line.find('GT') + 3])
        #         snp_count += 1
        #         fw.write(line)
        for line in fr:
            if line.startswith('##'):
                fw.write(line)
            else:
                break
        fw.write(line)

        # ii = 0
        vcf_reader = VCFReader(fpath=annotated_fpath)
        for line, record in tqdm(zip(fr, vcf_reader.get_record()), total=no_of_variants):
            if _SNP_condition(mode, SNP_type, record, effects):
                snp_count += 1
                fw.write(line)

            # ii += 1
            # if ii > 20:
            #     break

    print('Extracting SNPs from the annotated file finished. ' +
          f'We have extracted {snp_count} samples of the type {SNP_type}.')

    # print(effects)

    return output_fpath


def _SNP_condition(mode, SNP_type, record, effects):
    if mode == 'train':
        # if 'CLNSIG' in record['INFO'] and 'CLNVC' in record['INFO']:
        #     print(record['INFO']['CLNSIG'])
        #     print(record['INFO']['CLNVC'])
        return 'CLNVC' in record['INFO'] and \
                    record['INFO']['CLNVC'] == 'single_nucleotide_variant' \
                and 'MC' in record['INFO'] and any(
                    [role in record['INFO']['MC']
                        for role in interesting_roles[SNP_type]]) \
                and 'CLNSIG' in record['INFO'] and any(
                    [effect == record['INFO']['CLNSIG']
                        for effect in interesting_effects])
    else:
        if 'ANN' not in record['INFO']:
            return False
        ann = record['INFO']['ANN'].replace('nonsense_mediated_decay', '')
        # for role in interesting_roles[SNP_type]:
        #     if role in ann:
        #         effects[role] += 1
        return any([role in ann for role in interesting_roles[SNP_type]])
    # if mode == 'train':
    #     return (any([role in line for role in interesting_roles[SNP_type]]) \
    #         or 'nonsense' in line) \
    #             and any([effect in line for effect in interesting_effects])
    #         # For the training data, we only consider the examples with provided
    # else:
    #     return any([role in line for role in interesting_roles[SNP_type]]) 


def _labels_condition(record):
    return 'CLNSIG' in record['INFO'] \
        and record['INFO']['CLNSIG'] in interesting_effects  # interesting for training


def extract_hgvs_ids(fpath):

    flen = file_info.file_length(fpath)
    metalen = file_info.no_info_lines(fpath)
    no_of_variants = flen - metalen

    hgvs_fpath = results_fpath(fpath, 'hgvs')

    if os.path.exists(hgvs_fpath):
        print('HGVS ids have already been extracted.')
        return hgvs_fpath

    print('Extracting HGVS ids.')

    # vcf_reader = vcf.Reader(filename=fpath)
    vcf_reader = VCFReader(fpath=fpath)
    with open(fpath, 'r') as fr, open(hgvs_fpath, 'w') as fh:
        fh.write('ID\n')
        # for record in tqdm(vcf_reader, total=no_of_variants):
        for record in tqdm(vcf_reader.get_record(), total=no_of_variants):
            # for alt in record.ALT:
            #     if not (str(record.REF) in allowed_alleles and str(alt) in allowed_alleles):
            #         break
            #     fh.write(
            #         f'{_hgvs_id(record.CHROM, record.POS, record.REF, alt)}\n')
            for alt in record["ALT"]:
                for ref in record["REF"]:
                    if not(ref in allowed_alleles and alt in allowed_alleles):
                        break
                    fh.write(
                        f'{_hgvs_id(record["CHROM"], record["POS"], ref, alt)}\n')

    print('Extracting HGVS ids finished.')

    return hgvs_fpath


def extract_labels(fpath, mode, classification):
    if mode != 'train':
        print('We are not training the model. ' +
              'We do not need to extract the labels.')
        return None

    flen = file_info.file_length(fpath)
    metalen = file_info.no_info_lines(fpath)
    no_of_variants = flen - metalen

    labels_fpath = results_fpath(fpath, f'labels_{classification}')

    if os.path.exists(labels_fpath):
        print('The labels have already been extracted.')
        return labels_fpath

    print('Extracting the labels.')

    # vcf_reader = vcf.Reader(filename=fpath)
    vcf_reader = VCFReader(fpath=fpath)
    with open(fpath, 'r') as fr, open(labels_fpath, 'w') as fl:
        fl.write('ID, binary_label, multiclass_label\n')
        for record in tqdm(vcf_reader.get_record(), total=no_of_variants):
            # if not _labels_condition(record):
            #     break
            # for alt in record.ALT:
            #     if not(record.REF in allowed_alleles and alt in allowed_alleles):
            #         break
            #     fl.write(
            #         f'{_hgvs_id(record.CHROM, record.POS, record.REF, alt)},' +
            #             f'{descr_to_label[record.INFO["CLNSIG"][0]]}\n')
            for alt in record["ALT"]:
                for ref in record["REF"]: 
                    if not(ref in allowed_alleles and alt in allowed_alleles):
                        break
                    fl.write(
                        f'{_hgvs_id(record["CHROM"], record["POS"], ref, alt)},' +
                            f'{descr_to_label_binary[record["INFO"]["CLNSIG"]]},' +
                            f'{descr_to_label_multiclass[record["INFO"]["CLNSIG"]]}\n')

    print('Extracting the labels finished.')

    return labels_fpath