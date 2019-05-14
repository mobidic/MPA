#!/usr/bin/env python3
#
# Copyright (C) 2019
#

__author__ = 'Mobidic'
__authors__ = ['Henri Pegeot','Kevin Yauy','Charles Van Goethem','David Baux']
__copyright__ = 'Copyright (C) 2019'
__license__ = 'Academic License Agreement'
__version__ = '1.0.0'
__email__ = 'h-pegeot@chu-montpellier.fr'
__status__ = 'prod'

################################################################################
#Ó
# IMPORT
#
################################################################################
import vcf        # read vcf => PyVCF :https://pyvcf.readthedocs.io/en/latest/
import sys        # system command
import os         # os command
import re         # regex
import argparse   # for options
import logging    # logging messages
import subprocess # launch subprocess
import collections

########################################################################
#
# FUNCTIONS
#
########################################################################
def getSoftwarePath(software, expected_folder):
    """
    @summary: Returns the path to the software from the expected_folder if it is present or from the PATH environment variable.
    @param software: [str] Name of the software.
    @param expected_folder: [str] The folder where the software it is supposed to be in mSINGS.
    @return: [str] The path of the software.
    """
    path = os.path.join(expected_folder, software)  # Expected path in mSINGS directory
    if not os.path.exists(path):
        path = wich(software)  # Search in PATH
        if path is None:
            raise Exception("The software {} cannot be found in environment.".format(software))
    return path

def wich(software):
    """
    @summary: Returns the path to the software from the PATH environment variable.
    @param software: [str] Name of the software.
    @return: [str/None] The path of the software or None if it is not found.
    """
    soft_path = None
    PATH = os.environ.get('PATH')
    for current_folder in reversed(PATH.split(os.pathsep)):  # Reverse PATh elements to kept only the first valid folder
        eval_path = os.path.join(current_folder, software)
        if os.path.exists(eval_path):
            soft_path = eval_path
    return soft_path

def check_annotation(vcf_infos):
    """
    @summary: Chek if vcf followed the guidelines for annotations (17 are mandatory see full documentation)
    @param vcf_infos: [vcf.reader.infos] One record of the VCF
    @return: [None]
    """
    vcf_keys = [
        'Func.refGene',
    	'ExonicFunc.refGene',
    	'FATHMM_pred',
    	'dbscSNV_ADA_SCORE',
    	'dbscSNV_RF_SCORE',
    	#'dpsi_zscore',
        'spliceai_filtered',
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
    	'CLNSIG'
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

    log.debug("scores impact : " + str(scores_impact))

    for score, impact in scores_impact.items():
        if(impact == "D"):
            deleterious += 1
            available += 1
        elif(impact != None):
            available += 1

    if available > 0:
        score_adjusted = float(deleterious)/float(available) * 10


    log.debug(">> Return: ")
    log.debug({
        "adjusted":score_adjusted,
        "available":available,
        "deleterious":deleterious
    })
	# Return meta score and available tools
    return {
        "adjusted":score_adjusted,
        "available":available,
        "deleterious":deleterious
    }

# TODO: modulate clinvar score
def is_clinvar_pathogenic(clinsig):
    """
    @summary: Define if clinvar annotation predict this variant as pathogenic
    @param clinsig: [str] The clinvar annotation provided by the vcf
    @return: [int/bool] Rank (1) if is pathogenic and no Benign; False in other cases
    """
    # possible clinsig
    clinsig_possible = {
        "Benign": 0,
        "Benign/Likely_benign": 1,
        "Likely_benign": 2,
        "Uncertain_significance": 3,
        "Likely_pathogenic": 4,
        "Pathogenic/Likely_pathogenic": 5,
        "Pathogenic": 6,
        "other": None,
        "Affects": None,
        "drug_response": None,
        "confers_sensitivity": None,
        "risk_factor": None,
        "association": None,
        "association_not_found": None,
        "protective": None,
        "not_provided": None,
        "conflicting_data_from_submitters": None,
    }
    # No clinsig available
    if clinsig == None:
        return False

    # Test if "Pathogenic" or "Benign" match on clinsig
    match_pathogenic = re.search("pathogenic", clinsig, re.IGNORECASE)
    match_benign = re.search("benign", clinsig, re.IGNORECASE)
    match_conflicting = re.search("conflicting", clinsig, re.IGNORECASE)

    # Determine if clinvar as no doubt about pathogenicity
    if(match_pathogenic and not match_benign and not match_conflicting):
        return 1
    else:
        return False


def is_splice_impact(splices_scores, is_indel, funcRefGene):
    """
    @summary: Predict splicing effect of the variant
    @param splices_scores: [dict] The dictionnary of splicing scores
    @param is_indel: [bool] Boolean to define if variants is indel or not
    @param funcRefGene: [str] Annotation provided by refGene about the biological function
    @return: [int/bool] Rank (3,4,5,6,7 or 8) if is splicing impact; False in other cases
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
    # TODO: Replace for splice AI
    # Zscore_splice = (splices_scores["Zscore"] != None and
    #     splices_scores["ADA"] == None and
    #     splices_scores["RF"] == None and
    #     float(splices_scores["Zscore"]) < -2
    # )
    if(splices_scores["spliceAI"] != None):
        spliceAI_split = splices_scores["spliceAI"].split("\\x3b")
        spliceAI_annot = dict()
        for annot in spliceAI_split:
            annot_split = annot.split("\\x3d")
            spliceAI_annot[annot_split[0]] = annot_split[1]
    spliceAI_score_high = (splices_scores["spliceAI"] != None and
        (float(spliceAI_annot["DS_AG"]) > 0.8 or
        float(spliceAI_annot["DS_AL"]) > 0.8 or
        float(spliceAI_annot["DS_DG"]) > 0.8 or
        float(spliceAI_annot["DS_DL"]) > 0.8 )
    	)
    spliceAI_score_moderate = (splices_scores["spliceAI"] != None and
        (float(spliceAI_annot["DS_AG"]) > 0.5 or
        float(spliceAI_annot["DS_AL"]) > 0.5 or
        float(spliceAI_annot["DS_DG"]) > 0.5 or
        float(spliceAI_annot["DS_DL"]) > 0.5 )
    )
    spliceAI_score_low = (splices_scores["spliceAI"] != None and
        (float(spliceAI_annot["DS_AG"]) > 0.2 or
        float(spliceAI_annot["DS_AL"]) > 0.2 or
        float(spliceAI_annot["DS_DG"]) > 0.2 or
        float(spliceAI_annot["DS_DL"]) > 0.2 )
    )

    # Home made prediction of splice impact
    match_splicing = re.search("splicing", funcRefGene, re.IGNORECASE)
    home_splice = (is_indel and match_splicing)

    # Determine if there is a splicing impact
    if(RF_splice):
        return 3
    elif(ADA_splice):
        return 4
    # TODO: Replace for splice AI
    # elif(Zscore_splice):
    elif(spliceAI_score_high):
        return 5
    elif(spliceAI_score_moderate):
        return 6
    elif(spliceAI_score_low):
        return 9
    elif(home_splice):
        return 7
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
    match_nonframeshift = re.search("nonframeshift", exonicFuncRefGene, re.IGNORECASE)

    if(match_frameshift and not match_nonframeshift):
        return 2
    else:
        return False

def is_missense_impact(exonicFuncRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank (9) if is missense impact; False in other cases
    """
    match_missense = re.search("nonsynonymous_SNV", exonicFuncRefGene, re.IGNORECASE)

    if(match_missense):
        return 8
    else:
        return False

def is_unknown_impact(exonicFuncRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank (10) if is unknown impact; False in other cases
    """
    match_unknown = re.search("unknown", exonicFuncRefGene, re.IGNORECASE)

    if(match_unknown):
        return 10
    else:
        return False

################################################################################
#
# PROCESS
#
################################################################################
def main(args, logger):
    """
    @summary: Launch annotation with MPA score on a vcf.
    @param args: [Namespace] The namespace extract from the script arguments.
    param log: [Logger] The logger of the script.
    """
    global log
    log = logger

    # TODO: improve this ! already existing on pyVCF
    _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
    info_MPA_adjusted = _Info("MPA_adjusted", ".", "String", "MPA_adjusted : normalize MPA missense score from 0 to 10", "MPA", "0.3")
    info_MPA_available = _Info("MPA_available", ".", "String", "MPA_available : number of missense tools annotation available for this variant", "MPA", "0.3")
    info_MPA_deleterious = _Info("MPA_deleterious", ".", "String", "MPA_deleterious : number of missense tools that annotate this variant pathogenic", "MPA", "0.3")
    info_MPA_final_score = _Info("MPA_final_score", ".", "String", "MPA_final_score : unique score that take into account curated database, biological assumptions, splicing predictions and the sum of various predictors for missense alterations. Annotations are made for exonic and splicing variants up to +300nt.", "MPA", "0.3")
    info_MPA_impact = _Info("MPA_impact", ".", "String", "MPA_impact : pathogenic predictions (clinvar_pathogenicity, splice_impact, stop and frameshift_impact)", "MPA", "0.3")
    info_MPA_ranking = _Info("MPA_ranking", ".", "String", "MPA_ranking : prioritize variants with ranks from 1 to 7", "MPA", "0.3")

    with open(args.input, 'r') as f:
        log.info("Read VCF")
        vcf_reader = vcf.Reader(f)
        # TODO: improve this
        vcf_reader.infos.update({'MPA_adjusted':info_MPA_adjusted})
        vcf_reader.infos.update({'MPA_available':info_MPA_available})
        vcf_reader.infos.update({'MPA_deleterious':info_MPA_deleterious})
        vcf_reader.infos.update({'MPA_final_score':info_MPA_final_score})
        vcf_reader.infos.update({'MPA_impact':info_MPA_impact})
        vcf_reader.infos.update({'MPA_ranking':info_MPA_ranking})
        vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
        log.info("Check vcf annotations")

        try:
            check_annotation(vcf_reader.infos)
        except SystemExit as e:
            log.error(str(e))
            return 1

        log.info("Read the each variants")
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
            # TODO: Replace for splice AI
            splices_scores = {
                "ADA": record.INFO['dbscSNV_ADA_SCORE'][0],
                "RF": record.INFO['dbscSNV_RF_SCORE'][0],
                # "Zscore":record.INFO['dpsi_zscore'][0],
                "spliceAI":record.INFO['spliceai_filtered'][0],
            }

            # MPA aggregate the information to predict some effects
            meta_impact = {
                "clinvar_pathogenicity": False,
                "stop_impact": False,
                "splice_impact": False,
                "frameshift_impact": False,
                "unknown_impact": False
            }

            # Calculate adjusted score for each variants
            adjusted_score = calculate_adjusted_score(impacts_scores)

            # Determine if variant is well annotated with clinvar as deleterious
            meta_impact["clinvar_pathogenicity"] = is_clinvar_pathogenic(record.INFO['CLNSIG'][0])

            # Determine the impact on splicing
            meta_impact["splice_impact"] = is_splice_impact(splices_scores, record.is_indel, record.INFO['Func.refGene'][0])

            # Determine the exonic impact
            match_exonic = re.search("exonic", record.INFO['Func.refGene'][0], re.IGNORECASE)
            if(match_exonic and record.INFO['ExonicFunc.refGene'][0] != None):
                # Determine the stop impact
                meta_impact["stop_impact"] = is_stop_impact(record.INFO['ExonicFunc.refGene'][0])

                # Determine the frameshift impact
                meta_impact["frameshift_impact"] = is_frameshift_impact(record.INFO['ExonicFunc.refGene'][0])

                # Determine the missense impact
                meta_impact["missense_impact"] = is_missense_impact(record.INFO['ExonicFunc.refGene'][0])

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
                        if(impact == "unknown_impact" or impact == "missense_impact"):
                            adjusted_score["final_score"] = adjusted_score["adjusted"]
                        elif(impact == "splice_impact" and meta_impact["splice_impact"] == 6):
                            adjusted_score["final_score"] = 8
                        elif(impact == "splice_impact" and meta_impact["splice_impact"] == 9):
                            adjusted_score["final_score"] = 6
                        else:
                            adjusted_score["final_score"] = 10

            # if not ranking default value 10
            if not rank:
                rank = 10
                record.INFO['MPA_impact'] = "NULL,"
                adjusted_score["final_score"] = adjusted_score["adjusted"]

            log.debug("Ranking : " + str(rank))

            # write vcf output
            record.INFO['MPA_impact'] = record.INFO['MPA_impact'][:-1]
            record.INFO['MPA_ranking'] = rank
            for sc in adjusted_score:
                record.INFO['MPA_' + sc] = adjusted_score[sc]

            vcf_writer.write_record(record)
        vcf_writer.close()
