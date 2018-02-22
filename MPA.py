#!/usr/bin/env python2.7
#
# Copyright (C) 2018
#

__author__ = 'Henri Pegeot and Kevin Yauy and Charles Van Goethem'
__copyright__ = 'Copyright (C) 2018'
__license__ = 'Academic License Agreement'
__version__ = '0.1.0'
__email__ = 'h-pegeot@chu-montpellier.fr'
__status__ = 'dev'

################################################################################
#
# IMPORT
#
################################################################################
import vcf		# read vcf => PyVCF :https://pyvcf.readthedocs.io/en/latest/
import sys		# system command
import re		# regex
import argparse	# for options
import logging  # logging messages

########################################################################
#
# FUNCTIONS
#
########################################################################

def check_annotation(vcf_infos):
    """
    @summary: Chek if vcf followed the guidelines for annotations (17 are mandatory see full documentation)
    @param vcf_infos: [vcf.reader.infos] One record of the VCF
    @return: [None]
    """
    vcf_keys = [
    	'ExonicFunc.refGene',
    	'FATHMM_pred',
    	'Func.refGene',
    	'dbscSNV_ADA_SCORE',
    	'dbscSNV_RF_SCORE',
    	'dpsi_zscore',
    	'SIFT_pred',
    	'Polyphen2_HDIV_pred',
    	'Polyphen2_HVAR_pred',
    	'LRT_pred',
    	'MutationTaster_pred',
    	'FATHMM_pred',
    	'PROVEAN_pred',
    	'fathmm-MKL_coding_pred',
    	'MetaSVM_pred',
    	'MetaLR_pred',
    	'CLINSIG'
    ]

    if(not set(vcf_keys).issubset(vcf_infos)):
        sys.exit('VCF not correctly annotated. See documentation and provide a well annotated vcf (annotation with annovar).')

    return None

def check_split_variants(record):
    """
    @summary: Chek if vcf followed the specifications (only one reference) and guidelines pre-processed vcf (split variants)
    @param record: [vcf.model._record] One record of the VCF
    @return: [None]
    """
    if (len(str(record.REF).split(',')) > 1):
        sys.exit('Multi references on vcf. It seems that your vcf not followed the specifications.')

    if (len(record.ALT) > 1):
        sys.exit('Multi allelic variant on vcf. See documentation and provide a well processed vcf (split variants).')

    return None

def calculate_adjusted_score(scores_impact):
    """
    @summary: Calculate the adjusted score impact from 10 annotation score
    @param scores_impact: [dict] The dictionnary of impact score for a variant.
    @return: [dict] The dictionnary with adjusted, available and deleterious scores
    """
    deleterious = 0
    available = 0
    score_adjusted = 0

    for score, impact in scores_impact.items():
        if(impact == "D"):
            deleterious += 1
            available += 1
        elif(impact != "."):
            available += 1

    if available > 0:
        score_adjusted = deleterious/available * 10

	# Return meta score and available tools
    return {
        "adjusted":score_adjusted,
        "available":available,
        "deleterious":deleterious
    }

def is_clinvar_pathogenic(clinsig):
    """
    @summary: Define if clinvar annotation predict this variant as pathogenic
    @param clinsig: [str] The clinvar annotation provided by the vcf
    @return: [int/bool] Rank (1) if is pathogenic and no Benign; False in other cases
    """
    # No clinsig available
    if clinsig == None:
        return False

    # Test if "Pathogenic" or "Benign" match on clinsig
    match_pathogenic = re.search("pathogenic", clinsig, re.IGNORECASE)
    match_benign = re.search("benign", clinsig, re.IGNORECASE)

    # Determine if clinvar as no doubt about pathogenicity
    if(match_pathogenic and not match_benign):
        return 1
    else:
        return False


def is_splice_impact(splices_scores, is_indel, funcRefGene):
    """
    @summary: Predict splicing effect of the variant
    @param splices_scores: [dict] The dictionnary of splicing scores
    @param is_indel: [bool] Boolean to define if variants is indel or not
    @param funcRefGene: [str] Annotation provided by refGene about the biological function
    @return: [int/bool] Rank (3,4,5 or 6) if is splicing impact; False in other cases
    """

    # If ADA predict splicing impact
    ADA_splice = (splices_scores["ADA"] != None and
        float(splices_scores["ADA"]) >= 0.6
    )

    # If RF predict splicing impact
    RF_splice = (splices_scores["RF"] != None and
        float(splices_scores["RF"]) >= 0.6
    )

    # If Zscore predict splicing impact but no ADA and RF annotation
    Zscore_splice = (splices_scores["Zscore"] != None and
        splices_scores["ADA"] == None and
        splices_scores["RF"] == None and
        float(splices_scores["Zscore"]) < -2
    )

    # Home made prediction of splice impact
    match_splicing = re.search("splicing", funcRefGene, re.IGNORECASE)
    home_splice = (is_indel and match_splicing)

    # Determine if there is a splicing impact
    if(ADA_splice):
        return 4
    elif(RF_splice):
        return 3
    elif(Zscore_splice):
        return 5
    elif(home_splice):
        return 6
    else:
        return False

def is_stop_impact(exonicFuncRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [bool] Rank (2) if is stop impact; False in other cases
    """
    match_stoploss = re.search("stoploss", exonicFuncRefGene, re.IGNORECASE)
    match_stopgain = re.search("stopgain", exonicFuncRefGene, re.IGNORECASE)

    if(match_stopgain or match_stoploss):
        return 2
    else:
        return False

def is_frameshift_impact(exonicFuncRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank (2) if is frameshift impact; False in other cases
    """
    match_frameshift = re.search("frameshift", exonicFuncRefGene, re.IGNORECASE)

    if(match_frameshift):
        return 2
    else:
        return False

def is_unknown_impact(exonicFuncRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank (8) if is frameshift impact; False in other cases
    """
    match_frameshift = re.search("unknown", exonicFuncRefGene, re.IGNORECASE)

    if(match_frameshift):
        return 8
    else:
        return False

################################################################################
#
# PROCESS
#
################################################################################
def process(args, log):
    """
    @summary: Launch annotation with MPA score on a vcf.
    @param args: [Namespace] The namespace extract from the script arguments.
    param log: [Logger] The logger of the script.
    """

    with open(args.input, 'r') as f:
        log.info("Read VCF")
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(open('test.vcf', 'w'), vcf_reader)
        log.info("Check vcf annotations")

        try:
            check_annotation(vcf_reader.infos)
        except SystemExit as e:
            log.error(str(e))
            return

        log.info("Read the vcf")
        for record in vcf_reader:
            try:
                check_split_variants(record)
            except SystemExit as e:
                log.error(str(record))
                log.error(str(e))
                continue
            log.debug(str(record))

            # Deleterious impact scores
            impacts_scores = {
                "SIFT" : record.INFO['SIFT_pred'][0],
                "HDIV" : record.INFO['Polyphen2_HDIV_pred'][0],
                "HVAR" : record.INFO['Polyphen2_HVAR_pred'][0],
                "LRT" : record.INFO['LRT_pred'][0],
                "MutationTaster" : record.INFO['MutationTaster_pred'][0],
                "FATHMM" : record.INFO['FATHMM_pred'][0],
                "PROVEAN" : record.INFO['PROVEAN_pred'][0],
                "MKL" : record.INFO['fathmm-MKL_coding_pred'][0],
                "SVM" : record.INFO['MetaSVM_pred'][0],
                "LR" : record.INFO['MetaLR_pred'][0]
            }

            # Splicing impact scores
            splices_scores = {
                "ADA": record.INFO['dbscSNV_ADA_SCORE'][0],
                "RF": record.INFO['dbscSNV_RF_SCORE'][0],
                "Zscore":record.INFO['dpsi_zscore'][0],
            }

            # MPA aggregate the information to predict some effects
            meta_impact = {
                "clinvar_pathogenicity": False,
                "splice_impact": False,
                "stop_impact": False,
                "frameshift_impact": False,
                "unknown_impact": False
            }

            # Calculate adjusted score for each variants
            adjusted_score = calculate_adjusted_score(impacts_scores)

            # Determine if variant is well annotated with clinvar as deleterious
            meta_impact["clinvar_pathogenicity"] = is_clinvar_pathogenic(record.INFO['CLINSIG'][0])

            # Determine the impact on splicing
            meta_impact["splice_impact"] = is_splice_impact(splices_scores, record.is_indel, record.INFO['Func.refGene'][0])

            # Determine the exonic impact
            match_exonic = re.search("exonic", record.INFO['Func.refGene'][0], re.IGNORECASE)
            if(match_exonic and record.INFO['ExonicFunc.refGene'][0] != None):
                # Determine the stop impact
                meta_impact["stop_impact"] = is_stop_impact(record.INFO['ExonicFunc.refGene'][0])

                # Determine the frameshift impact
                meta_impact["frameshift_impact"] = is_frameshift_impact(record.INFO['ExonicFunc.refGene'][0])

                # Determine if unknown impact (misunderstand gene)
                # NOTE: /!\ Be careful to updates regularly your databases /!\
                meta_impact["unknown_impact"] = is_unknown_impact(record.INFO['ExonicFunc.refGene'][0])
            log.debug("Meta score : " + str(meta_impact))

            # Ranking of variants
            rank = False
            record.INFO['MPA_impact'] = ""
            for impact in meta_impact:
                if (meta_impact[impact]):
                    record.INFO['MPA_impact'] = record.INFO['MPA_impact'] + impact + ","
                    if(meta_impact[impact]<rank or not rank):
                        rank = meta_impact[impact]
                        if(impact == "unknown_impact"):
                            adjusted_score["final_score"] = adjusted_score["adjusted"]
                        else:
                            adjusted_score["final_score"] = 10
            if not rank:
                rank = 7
                record.INFO['MPA_impact'] = "NULL,"
                adjusted_score["final_score"] = adjusted_score["adjusted"]

            log.debug("Ranking : " + str(rank))

            record.INFO['MPA_impact'] = record.INFO['MPA_impact'][:-1]
            record.INFO['MPA_ranking'] = rank
            for sc in adjusted_score:
                record.INFO['MPA_' + sc] = adjusted_score[sc]

            vcf_writer.write_record(record)




################################################################################
#
# CLASS
#
################################################################################
class LoggerAction(argparse.Action):
    """
    @summary: Manages logger level parameters (The value "INFO" becomes logging.info and so on).
    """
    def __call__(self, parser, namespace, values, option_string=None):
        log_level = None
        if values == "DEBUG":
            log_level = logging.DEBUG
        elif values == "INFO":
            log_level = logging.INFO
        elif values == "WARNING":
            log_level = logging.WARNING
        elif values == "ERROR":
            log_level = logging.ERROR
        elif values == "CRITICAL":
            log_level = logging.CRITICAL
        setattr(namespace, self.dest, log_level)

################################################################################
#
# MAIN
#
################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Annotate VCF with Mobidic Prioritization Algorithm score (MPA).")
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)

    group_input = parser.add_argument_group('Inputs') # Inputs
    group_input.add_argument('-i', '--input', required=True, help="The vcf file to annotate (format: VCF). This vcf must be annotate with annovar.")

    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output', required=True, help="The output vcf file with annotation (format : VCF)")
    args = parser.parse_args()

    # Process
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger("MPA_score")
    log.setLevel(args.logging_level)
    log.info("Start MPA annotation")
    log.info("Command: " + " ".join(sys.argv))
    process(args, log)
    log.info("End MPA annotation")
