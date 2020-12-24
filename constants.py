META_DATA = '##'
HEADER = '#'


NS_SNP = 'nsSNP'
S_SNP = 'sSNP'
I_SNP = 'iSNP'


CLNSIG = {
    0: 'Uncertain_significance',
    1: 'Uncertain_significance',
    2: 'Uncertain_significance',
    3: 'Uncertain_significance',
    4: 'Uncertain_significance'
}


fromatdesc2code = {
    'Read depth for each allele': 'AD',
    'Read depth for each allele on the forward strand': 'ADF',
    'Read depth for each allele on the reverse strand': 'ADR',
    'Read depth': 'DP',
    'Expected alternate allele counts': 'EC',
    'Filter indicating if this genotype was “called”': 'FT',
    'Genotype likelihoods': 'GL',
    # 'Conditional genotype quality': 'GP',
    'Conditional genotype quality': 'GQ',
    'Genotype': 'GT',
    'Haplotype quality': 'HQ',
    'RMS mapping quality': 'MQ',
    'Phred-scaled genotype likelihoods rounded to the closest integer': 'PL',
    'Phasing quality': 'PQ',
    'Phase set': 'PS'
}


interesting_effects = {
    'Benign',
    'Benign/Likely_benign',
    'Likely_benign',
    'Likely_pathogenic',
    'Pathogenic/Likely_pathogenic',
    'Pathogenic'
}


interesting_roles = {
    'nsSNP': {
        'missense_variant',
        'start_lost',
        'stop_lost',
        'start_gained',
        'stop_gained',
        'nonsense',
        'splice_region_variant'
    },
    'sSNP': {
        'synonymous_variant',
    },
    'iSNP': {
        'intron_variant',
    }
}


all_roles_snpEff = {
    '3_prime_UTR_variant',
    '5_prime_UTR_premature_start_codon_gain_variant',
    '5_prime_UTR_variant',
    'downstream_gene_variant',
    'initiator_codon_variant',
    'initiator_codon_variant&splice_region_variant',
    'intergenic_region',
    'intragenic_variant',
    'intron_variant',
    'missense_variant',
    'missense_variant&splice_region_variant',
    'non_coding_transcript_exon_variant',
    'splice_acceptor_variant&intron_variant',
    'splice_acceptor_variant&splice_donor_variant&intron_variant',
    'splice_acceptor_variant&splice_region_variant&intron_variant',
    'splice_donor_variant&intron_variant',
    'splice_region_variant',
    'splice_region_variant&intron_variant',
    'splice_region_variant&non_coding_transcript_exon_variant',
    'splice_region_variant&stop_retained_variant',
    'splice_region_variant&synonymous_variant',
    'start_lost',
    'start_lost&splice_region_variant',
    'stop_gained',
    'stop_gained&splice_region_variant',
    'stop_lost',
    'stop_lost&splice_region_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'upstream_gene_variant'
}


allowed_alleles = {'A', 'T', 'C', 'G'}


descr_to_label_binary = {
    'Benign': 0,
    'Benign/Likely_benign': 0,
    'Likely_benign': 0,
    'Likely_pathogenic': 1,
    'Pathogenic/Likely_pathogenic': 1,
    'Pathogenic': 1
}


descr_to_label_multiclass = {
    'Benign': 0,
    'Benign/Likely_benign': 1,
    'Likely_benign': 1,
    'Likely_pathogenic': 2,
    'Pathogenic/Likely_pathogenic': 2,
    'Pathogenic': 3
}


label_to_descr_binary = {
    0: 'Benign',
    1: 'Pathogenic'
}


label_to_descr_multiclass = {
    0: 'Benign',
    1: 'Likely_benign',
    2: 'Likely_pathogenic',
    3: 'Pathogenic'
}


n2c = {
    'A': 0,
    'T': 1,
    'C': 2,
    'G': 3,
}


default_mv_fields = {
    'nsSNP': [
        'cadd.gc',
        'cadd.phylop',
        'cadd.phast_cons',
        'cadd.dst2splice',
        'cadd.min_dist_tss',
        'cadd.min_dist_tse',
        'cadd.encode.sig',
        'cadd.encode',
        'cadd.gerp',
        'dbnsfp.uniprot.acc',
        'snpeff.ann.hgvs_p'
    ],
    'sSNP': [
        'cadd.phylop',
        'cadd.phast_cons',
        'cadd.gc',
        'cadd.cpg',
        'cadd.min_dist_tse',
        'cadd.min_dist_tss',
        'cadd.gerp',
        'cadd.encode',
        'dbnsfp.uniprot.acc',
        'snpeff.ann.hgvs_p'
    ],
    'iSNP': [
        'cadd.phylop',
        'cadd.phast_cons',
    ]
}

