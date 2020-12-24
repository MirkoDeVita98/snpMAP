from tqdm import tqdm

from utils import parser, file_info
from constants import META_DATA, HEADER, fromatdesc2code, label_to_descr_multiclass, label_to_descr_binary
from preprocessing.preprocess import results_fpath, _hgvs_id
from io_utils.reading import VCFReader


REF_COLUMN = 3
ALT_COLUMN = 4
INFO_COLUMN = 7

config = parser.load_config()
results_dir = config['output']['results_dir']
classification = config['classification']


class ResultsWriter:


    def __init__(self, fpath):
        self.input_fpath = fpath
        self.results_fpath = results_fpath(
            fpath, 'results', ending='vcf', folder=results_dir)
        self.file = None
        self.columns = None
        self.columns_full = None
        self.formats = []
        self.infos = []
        self.flag_infos = set()


    def write_results(self, results):

        print('Writing out the results.')

        conv = label_to_descr_binary if classification == 'binary' else label_to_descr_multiclass

        flen = file_info.file_length(self.input_fpath)
        metalen = file_info.no_info_lines(self.input_fpath)
        no_of_variants = flen - metalen

        ids = set(results.keys())

        with open(self.input_fpath, 'r') as fr, \
            open(self.results_fpath, 'w') as fw:
            for line in fr:
                if line.startswith('##'):
                    fw.write(line)
                else:
                    break
            fw.write(
                '##INFO=<ID=PREDICTION,Number=.,Type=String,Description="Classification of the variant by out model">\n')
            fw.write(line)

            vcf_reader = VCFReader(fpath=self.input_fpath)
            for line, record in tqdm(
                zip(fr, vcf_reader.get_record()), total=no_of_variants):
                for alt in record["ALT"]:
                    for ref in record["REF"]:
                        line_parts = line.split('\t')
                        hgvs = _hgvs_id(record["CHROM"], record["POS"], ref, alt)
                        if hgvs in ids:
                            line_parts[INFO_COLUMN] = f'PREDICTION={conv[results[hgvs]]};{line_parts[INFO_COLUMN]}'
                        else:
                            line_parts[INFO_COLUMN] = f'PREDICTION=Unknown;{line_parts[INFO_COLUMN]}'
                        fw.write('\t'.join(line_parts))
        
        print('Writing out the results finished.')
