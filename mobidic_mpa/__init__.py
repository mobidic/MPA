#!/usr/bin/env python3
#
# Copyright (C) 2017-2022
#

__author__ = 'Mobidic'
__authors__ = [
    'Henri Pegeot',
    'Kevin Yauy',
    'Charles Van Goethem',
    'Thomas Guignard',
    'David Baux'
]
__copyright__ = 'Copyright (C) 2017-2022'
__license__ = 'Academic License Agreement'
__version__ = '1.2.4'
__email__ = 'c-vangoethem@chu-montpellier.fr'
__status__ = 'prod'

###############################################################################
#
# IMPORT
#
###############################################################################
import vcf        # read vcf => PyVCF :https://pyvcf.readthedocs.io/en/latest/
import sys        # system command
import re         # regex
import collections
import tqdm


########################################################################
#
# FUNCTIONS
#
########################################################################
def check_annotation(vcf_infos, no_refseq_version=True):
    """
    @summary: Chek if vcf followed the guidelines for annotations (17 are \
        mandatory see full documentation)
    @param vcf_infos: [vcf.reader.infos] One record of the VCF
    @return: [None]
    """

    refSeqExt = 'refGene' if no_refseq_version else 'refGeneWithVer'
    vcf_keys = [
        'Func.{}'.format(refSeqExt),
        'ExonicFunc.{}'.format(refSeqExt),
        'dbscSNV_ADA_SCORE',
        'dbscSNV_RF_SCORE',
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

    log.debug(vcf_keys)

    if(not set(vcf_keys).issubset(vcf_infos)):
        sys.exit(
            "VCF not correctly annotated. See documentation and provide "
            "a well annotated vcf (annotation with annovar)."
        )

    return None


def check_split_variants(record):
    """
    @summary: Chek if vcf followed the specifications (only one reference) \
        and guidelines pre-processed vcf (split variants)
    @param record: [vcf.model._record] One record of the VCF
    @return: [None]
    """
    if (len(str(record.REF).split(',')) > 1):
        sys.exit(
            "Multi references on vcf. It seems that your vcf not "
            "followed the specifications."
        )

    if (len(record.ALT) > 1):
        sys.exit(
            "Multi allelic variant on vcf. See documentation and provide "
            "a well processed vcf (split variants)."
        )

    return None


def calculate_adjusted_score(scores_impact):
    """
    @summary: Calculate the adjusted score impact from 10 annotation score
    @param scores_impact: [dict] The dictionnary of impact score for a variant.
    @return: [dict] The dictionnary with adjusted, available and deleterious \
        scores
    """
    deleterious = 0
    available = 0
    score_adjusted = 0

    log.debug(f"scores impact : {scores_impact}")

    for score, impact in scores_impact.items():
        if(impact == "D" or impact == "A"):
            deleterious += 1
            available += 1
        elif(impact is not None):
            available += 1

    if available > 0:
        score_adjusted = float(deleterious)/float(available) * 10

    log.debug(">> Return: ")
    log.debug({
        "adjusted": score_adjusted,
        "available": available,
        "deleterious": deleterious
    })

    # Return meta score and available tools
    return {
        "adjusted": score_adjusted,
        "available": available,
        "deleterious": deleterious
    }


# TODO: modulate clinvar score
def is_clinvar_pathogenic(clinsig):
    """
    @summary: Define if clinvar annotation predict this variant as pathogenic
    @param clinsig: [str] The clinvar annotation provided by the vcf
    @return: [int/bool] Rank (1) if is pathogenic and no Benign; False in \
        other cases
    """
    # No clinsig available
    if clinsig is None:
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
    @param funcRefGene: [str] Annotation provided by refGene about the \
        biological function
    @return: [int/bool] Rank (3,4,5,6,7 or 8) if is splicing impact; False in \
        other cases
    """

    # If ADA predict splicing impact
    ADA_splice = (
        splices_scores["ADA"] is not None and
        float(splices_scores["ADA"]) >= 0.6
    )

    # If RF predict splicing impact
    RF_splice = (
        splices_scores["RF"] is not None and
        float(splices_scores["RF"]) >= 0.6
    )

    # If Zscore predict splicing impact but no ADA and RF annotation
    if(splices_scores["spliceAI"] is not None):
        spliceAI_split = splices_scores["spliceAI"].split("\\x3b")
        spliceAI_annot = dict()
        for annot in spliceAI_split:
            annot_split = annot.split("\\x3d")
            if len(annot_split) > 1:
                spliceAI_annot[annot_split[0]] = annot_split[1]

    spliceAI_score_high = (
        splices_scores["spliceAI"] is not None and
        (
            float(spliceAI_annot["DS_AG"]) > 0.8 or
            float(spliceAI_annot["DS_AL"]) > 0.8 or
            float(spliceAI_annot["DS_DG"]) > 0.8 or
            float(spliceAI_annot["DS_DL"]) > 0.8
        )
    )
    spliceAI_score_moderate = (
        splices_scores["spliceAI"] is not None and
        (
            float(spliceAI_annot["DS_AG"]) > 0.5 or
            float(spliceAI_annot["DS_AL"]) > 0.5 or
            float(spliceAI_annot["DS_DG"]) > 0.5 or
            float(spliceAI_annot["DS_DL"]) > 0.5
        )
    )
    spliceAI_score_low = (
        splices_scores["spliceAI"] is not None and
        (
            float(spliceAI_annot["DS_AG"]) > 0.2 or
            float(spliceAI_annot["DS_AL"]) > 0.2 or
            float(spliceAI_annot["DS_DG"]) > 0.2 or
            float(spliceAI_annot["DS_DL"]) > 0.2
        )
    )

    # Home made prediction of splice impact
    match_splicing = re.search("splicing", funcRefGene, re.IGNORECASE)
    home_splice = (is_indel and match_splicing)

    # Determine if there is a splicing impact
    if(RF_splice):
        return 3
    elif(ADA_splice):
        return 3
    elif(spliceAI_score_high):
        return 4
    elif(spliceAI_score_moderate):
        return 6
    elif(spliceAI_score_low):
        return 8
    elif(home_splice):
        return 8
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


def is_indel_impact(exonicFuncRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank (2) if is frameshift impact; False in other cases
    """
    match_frameshift = re.search(
        "frameshift",
        exonicFuncRefGene,
        re.IGNORECASE)
    match_nonframeshift = re.search(
        "nonframeshift",
        exonicFuncRefGene,
        re.IGNORECASE)

    if(match_frameshift and not match_nonframeshift):
        return 2
    elif(match_nonframeshift):
        return 8
    else:
        return False


def is_missense_impact(exonicFuncRefGene, adjusted_score):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank () if is missense impact; False in other cases
    """
    match_missense = re.search(
        "nonsynonymous_SNV",
        exonicFuncRefGene,
        re.IGNORECASE)

    if(match_missense):
        if(adjusted_score > 6):
            return 5
        elif(adjusted_score > 2):
            return 7
        else:
            return 9
    else:
        return False


def is_unknown_impact(exonicFuncRefGene):
    """
    @summary: if no effect known
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [int/bool] Rank (10) if is unknown impact; False in other cases
    """
    match_unknown = re.search("unknown", exonicFuncRefGene, re.IGNORECASE)

    if(match_unknown):
        return 10
    else:
        return False


###############################################################################
#
# PROCESS
#
###############################################################################
def main(args, logger):
    """
    @summary: Launch annotation with MPA score on a vcf.
    @param args: [Namespace] The namespace extract from the script arguments.
    param log: [Logger] The logger of the script.
    """
    global log
    log = logger

    # TODO: improve this ! already existing on pyVCF
    _Info = collections.namedtuple(
        'Info',
        ['id', 'num', 'type', 'desc', 'source', 'version']
    )
    info_MPA_adjusted = _Info(
        "MPA_adjusted",
        ".",
        "String",
        "MPA_adjusted : normalize MPA missense score from 0 to 10",
        "MPA",
        __version__
    )
    info_MPA_available = _Info(
        "MPA_available",
        ".",
        "String",
        "MPA_available : number of missense tools annotation available for this variant",
        "MPA",
        __version__
    )
    info_MPA_deleterious = _Info(
        "MPA_deleterious",
        ".",
        "String",
        "MPA_deleterious : number of missense tools that annotate this variant pathogenic",
        "MPA",
        __version__
    )
    info_MPA_final_score = _Info(
        "MPA_final_score",
        ".", "String",
        "MPA_final_score : unique score that take into account curated \
        database, biological assumptions, splicing predictions and the sum of \
        various predictors for missense alterations. Annotations are made for \
        exonic and splicing variants up to +300nt.",
        "MPA",
        __version__
    )
    info_MPA_impact = _Info(
        "MPA_impact",
        ".",
        "String",
        "MPA_impact : pathogenic predictions (clinvar_pathogenicity, splice_impact, stop, frameshift_impact & indel_impact)",
        "MPA",
        __version__
    )
    info_MPA_ranking = _Info(
        "MPA_ranking",
        ".",
        "String",
        "MPA_ranking : prioritize variants with ranks from 1 to 10",
        "MPA",
        __version__
    )
    refSeqExt = 'refGene' if args.no_refseq_version else 'refGeneWithVer'

    with open(args.input, 'r') as f:
        log.info("Read VCF file")
        vcf_reader = vcf.Reader(f)
        count = -1
        if not args.no_progress_bar:
            count = sum(1 for _ in vcf_reader)
            log.info(f"Number of variants : {count}")
            f.seek(0)
            vcf_reader = vcf.Reader(f)

        # TODO: improve this
        vcf_reader.infos.update({'MPA_adjusted': info_MPA_adjusted})
        vcf_reader.infos.update({'MPA_available': info_MPA_available})
        vcf_reader.infos.update({'MPA_deleterious': info_MPA_deleterious})
        vcf_reader.infos.update({'MPA_final_score': info_MPA_final_score})
        vcf_reader.infos.update({'MPA_impact': info_MPA_impact})
        vcf_reader.infos.update({'MPA_ranking': info_MPA_ranking})
        vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
        log.info("Check vcf annotations")

        try:
            check_annotation(vcf_reader.infos, args.no_refseq_version)
        except SystemExit as e:
            log.error(str(e))
            sys.exit(1)

        log.info("Read each variants")
        for record in tqdm.tqdm(vcf_reader, total=count):
            log.debug(str(record))

            try:
                check_split_variants(record)
            except SystemExit as e:
                log.error(str(record))
                log.error(str(e))
                sys.exit(2)

            # Deleterious impact scores
            impacts_scores = {
                "SIFT": record.INFO['SIFT_pred'][0],
                "HDIV": record.INFO['Polyphen2_HDIV_pred'][0],
                "HVAR": record.INFO['Polyphen2_HVAR_pred'][0],
                "LRT": record.INFO['LRT_pred'][0],
                "MutationTaster": record.INFO['MutationTaster_pred'][0],
                "FATHMM": record.INFO['FATHMM_pred'][0],
                "PROVEAN": record.INFO['PROVEAN_pred'][0],
                "MKL": record.INFO['fathmm-MKL_coding_pred'][0],
                "SVM": record.INFO['MetaSVM_pred'][0],
                "LR": record.INFO['MetaLR_pred'][0]
            }

            # Splicing impact scores
            splices_scores = {
                "ADA": record.INFO['dbscSNV_ADA_SCORE'][0],
                "RF": record.INFO['dbscSNV_RF_SCORE'][0],
                "spliceAI": record.INFO['spliceai_filtered'][0],
            }

            # MPA aggregate the information to predict some effects
            meta_impact = {
                "clinvar_pathogenicity": False,
                "stop_impact": False,
                "splice_impact": False,
                "frameshift_impact": False,
                "indel_impact": False,
                "unknown_impact": False
            }

            # Calculate adjusted score for each variants
            adjusted_score = calculate_adjusted_score(impacts_scores)

            # Determine if variant is annotated with clinvar as deleterious
            meta_impact["clinvar_pathogenicity"] = is_clinvar_pathogenic(
                record.INFO['CLNSIG'][0]
            )

            FuncKey = f'Func.{refSeqExt}'
            ExonicFuncKey = f'ExonicFunc.{refSeqExt}'

            # Determine the impact on splicing
            meta_impact["splice_impact"] = is_splice_impact(
                splices_scores,
                record.is_indel,
                record.INFO[FuncKey][0]
            )

            # Determine the exonic impact
            match_exonic = re.search(
                "exonic",
                record.INFO[FuncKey][0],
                re.IGNORECASE
            )
            if (
                match_exonic and
                record.INFO[ExonicFuncKey][0] is not None
            ):
                # Determine the stop impact
                meta_impact["stop_impact"] = is_stop_impact(
                    record.INFO[ExonicFuncKey][0])

                # Determine the frameshift impact
                if is_indel_impact(record.INFO[ExonicFuncKey][0]) == 8:
                    meta_impact["indel_impact"] = 8
                if is_indel_impact(record.INFO[ExonicFuncKey][0]) == 2:
                    meta_impact["frameshift_impact"] = 2

                # Determine the missense impact
                meta_impact["missense_impact"] = is_missense_impact(
                    record.INFO[ExonicFuncKey][0],
                    adjusted_score["adjusted"])

                # Determine if unknown impact (misunderstand gene)
                # NOTE: /!\ Be careful to updates regularly your databases /!\
                meta_impact["unknown_impact"] = is_unknown_impact(
                    record.INFO[ExonicFuncKey][0])

            log.debug(f"Meta score : {meta_impact}")

            # Ranking of variants
            rank = False
            record.INFO['MPA_impact'] = ""
            for impact in meta_impact:
                if (meta_impact[impact]):
                    record.INFO['MPA_impact'] = (
                        f"{record.INFO['MPA_impact']}"
                        f"{impact},"
                    )

                    if(meta_impact[impact] < rank or not rank):

                        rank = meta_impact[impact]

                        if (
                            impact == "unknown_impact" or
                            impact == "missense_impact"
                        ):
                            adjusted_score["final_score"] = \
                                adjusted_score["adjusted"]
                        elif (
                            impact == "splice_impact" and
                            meta_impact["splice_impact"] == 6
                        ):
                            adjusted_score["final_score"] = 6
                        elif (
                            impact == "splice_impact" and
                            meta_impact["splice_impact"] == 8
                        ):
                            adjusted_score["final_score"] = 2
                        elif (
                            impact == "indel_impact" and
                            meta_impact["indel_impact"] == 8
                        ):
                            adjusted_score["final_score"] = 8
                        elif (
                            impact == "frameshift_impact" and
                            meta_impact["frameshift_impact"] == 2
                        ):
                            adjusted_score["final_score"] = 2
                        else:
                            adjusted_score["final_score"] = 10

            # if not ranking default value 10
            if not rank:
                rank = 10
                record.INFO['MPA_impact'] = "NULL,"
                adjusted_score["final_score"] = adjusted_score["adjusted"]

            log.debug(f"Ranking : {rank}")

            # write vcf output
            record.INFO['MPA_impact'] = record.INFO['MPA_impact'][:-1]
            record.INFO['MPA_ranking'] = rank
            for sc in adjusted_score:
                record.INFO['MPA_' + sc] = adjusted_score[sc]

            vcf_writer.write_record(record)
        vcf_writer.close()
