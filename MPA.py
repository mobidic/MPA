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


import vcf		# read vcf => PyVCF :https://pyvcf.readthedocs.io/en/latest/
import sys		# system command
import csv		# read, write csv
import re		# regex
import argparse	# for options
import os		# for options
import pprint	# pretty print
import logging

# ==============================================================================

########################################################################
#
# FUNCTIONS
#
########################################################################

def score(exonicFunction = ".", ref = ".", alt = ".", function = ".",
		ADAscore = ".", RFscore = ".", Zscore = ".", SIFT = ".", HDIV = ".",
		HVAR = ".", LRT = ".", MutationTaster = ".", FATHMM = ".",
		PROVEAN = ".", MKL = ".", SVM = ".", LR = ".", CLINSIG = ".",):

	"""
	Calculate a score based on pathogenicity

    Keyword arguments:
	    exonicFunction	-- ExonicFunc_refGene (default ".")
		ref				-- ref (default ".")
		alt				-- alt (default ".")
		function		-- Func_refGene (default ".")
		ADAscore		-- SNV_ADA_SCORE (default ".")
		RFscore			-- dbscSNV_RF_SCORE (default ".")
		Zscore			-- dpsi_zscore (default ".")
		SIFT			-- SIFT_pred (default ".")
		HDIV			-- Polyphen2_HDIV_pred (default ".")
		HVAR			-- Polyphen2_HVAR_pred (default ".")
		LRT				-- LRT_pred (default ".")
		MutationTaster	-- MutationTaster_pred (default ".")
		FATHMM			-- FATHMM_pred (default ".")
		PROVEAN			-- PROVEAN_pred (default ".")
		MKL				-- fathmm_MKL_coding_pred (default ".")
		SVM				-- MetaSVM_pred (default ".")
		LR				-- MetaLR_pred (default ".")

	Return tuple 'score','level' (score max = 10, min = 0):
		'10db','clinvar'	: deleterious impact (10) reported as so in ClinVar
		'10sfs','stop'		: deleterious impact (10) by stop gain or loss
		'10sfs','frameshift': deleterious impact (10) by frameshift indel
		'10sp','splice'		: deleterious impact (10) on splicing
		'u','na'			: unknown exonic function
		'na','0'			: cannot predict impact ; no scores prediction
		'3.3','3/9'			: score of impact (3.3) ; 3 deleterious score on 9
    """

	# Create a dictionnary of impact score tools used
	scores_impact = {"SIFT" : SIFT, "HDIV" : HDIV, "HVAR" : HVAR, "LRT" : LRT,
		"MutationTaster" : MutationTaster, "FATHMM" : FATHMM,
		"PROVEAN" : PROVEAN, "MKL" : MKL, "SVM" : SVM, "LR" : LR}
# for record in vcf_reader:
# 	diff = ((len(record.REF)) - (len(str(record.ALT[0]))))
#
# 	if diff %3 != 0 and record.INFO['Func.refGene'][0] == "exonic" :
# 		print ("diff : ", diff)
# 		print ("Func.refGene : ", record.INFO['Func.refGene'][0])
# 		print ("ExonicFunc.refGene : ", record.INFO['ExonicFunc.refGene']),
# 		print ("longueur REF :", len(record.REF), "longueur ALT :", len(str((record.ALT[0]))))
# 		print ("(lg.REF - lg.ALT) :", diff)
# 		print ("REF : ", record.REF, "ALT : ", record.ALT)
# 		print (str(record.INFO['Gene.refGene'][0]))
	# Test and calculate the score.
	deleterious = 0
	available = 0
	score_adjusted = 0

	for score, impact in scores_impact.items():
		if(impact == "D"):
			deleterious += 1
			available += 1
		elif(impact != "."):
			available += 1

	# Return meta score and available tools
	if available > 0:
		score_adjusted = deleterious/available * 10
	# Select variant described by ClinVar database entries as "P/pathogenic"
	# and not as "B/benign"
	if ("Pathogenic" in str(CLINSIG) or "pathogenic" in str(CLINSIG)) and \
	((not "Benign" in str(CLINSIG)) and (not "benign" in str(CLINSIG))):
		return(score_adjusted, "10db", available, "clinvar",1)
	# Test impact on stop codon

	if (exonicFunction =="stopgain" or exonicFunction =="stoploss"
	or exonicFunction == "frameshift_deletion" or exonicFunction =="frameshift_insertion"):
		return(score_adjusted, "10sfs", available, "stop",2)

	if ((len(str(ref))) - (len(str(alt)))) % 3 !=0 and exonicFunction == "exonic":
		return(score_adjusted, "10sfs", available, "frameshift",2)

	# Test impact on splice function
	if(str(ADAscore) != "."):
		if(float(ADAscore) >= 0.6):
			return(score_adjusted, "10spADA",available, "splice",4)
	if(str(RFscore) != "."):
		if(float(RFscore) >= 0.6):
			return(score_adjusted, "10spRF",available, "splice",3)
	if(str(Zscore) != "."):
		if(float(Zscore) < -2):
			if (ADAscore == ".") and (RFscore == ".") :
				return(score_adjusted, "10sp", available, "splice",5)
	if (((ref == "-") or (alt == "-")) and re.match("splicing",function)):
		return(score_adjusted, "10sp_indel",available, "splice",6)


	# Test if variant maps to multiple location, return "u" (UNKNOWN) if it does
	if(exonicFunction == "unknown"):
		return(score_adjusted, "u",available, "u",8)
	#return 'na' if no tools score are available
	#return deleterious/available * 10, "%s/%s" % (deleterious, available)
	else:
		return (score_adjusted, score_adjusted, available, "na",7)

# # ==============================================================================
# # ==============================================================================
#
# # Writing the header
# head = "Rank\tScore\tGene.refGene\t"
#
# # Writing sample name in the header
# for elt in vcf_reader.samples:
# 	head = head + elt + "\t"
# # Writing the header (second part)
#
#
# head = head + "ExAc\tClinSig\tFunc.refGene\t\
# ExonicFunc.refGene\tAAChange.refGene\t#CHROM\tPOS\tREF\tALT\tADA\tRF\t\
# Spidex\tScore\tNumber of tools"
#
# # ==============================================================================
#
# # Iterating through the input VCF
# for record in vcf_reader:
#
# # =======																========
# # Testing if the input vcf is correctly annotated
#
# 	error = False
# 	try:
# 		record.REF 	= ["." if record.REF is None else record.REF]
# 	except KeyError:
# 		print("The input VCF file provided lacks the annotation by REF.")
# 		error=True
#
# 	try:
# 		record.ALT 	= ["." if elt is None else elt for elt in record.ALT]
# 	except KeyError:
# 		print("The input VCF file provided lacks the annotation by ALT.")
# 		error=True
#
# 	for elt in vcf_key:
# 		error = test_vcf_annotation(record, elt, error)
#
# # 	if error==True :
# 		sys.exit("\n" + "Please refer yourself to the documentation in order to\
#  process VCF in a valid format" + "\n")
#
# # =======                                                              =========
#
# # call score function, instanciate it to "s"
#
# 	s = score(
# 		record.INFO['ExonicFunc.refGene'][0],
# 		record.REF,record.ALT[0],
# 		record.INFO['Func.refGene'][0],
# 		record.INFO['dbscSNV_ADA_SCORE'][0],
# 		record.INFO['dbscSNV_RF_SCORE'][0],
# 		record.INFO['dpsi_zscore'][0],
# 		record.INFO['SIFT_pred'][0],
# 		record.INFO['Polyphen2_HDIV_pred'][0],
# 		record.INFO['Polyphen2_HVAR_pred'][0],
# 		record.INFO['LRT_pred'][0],
# 		record.INFO['MutationTaster_pred'][0],
# 		record.INFO['FATHMM_pred'][0],
# 		record.INFO['PROVEAN_pred'][0],
# 		record.INFO['fathmm-MKL_coding_pred'][0],
# 		record.INFO['MetaSVM_pred'][0],
# 		record.INFO['MetaLR_pred'][0],
# 		record.INFO['CLINSIG']
# 	)
#
# # ==============================================================================
#
# # Writing output csv
# 	geno	= ""
# 	for sample in record.samples:
# 			geno += ("'" + sample['GT'] + "\t")
# 	f.write(
# 		str(s[4]) + "\t" +
# 		str(s[1]) + "\t" +
# 		str(record.INFO['Gene.refGene'][0]) +"\t" +
# 		str(geno) +
# 		str(record.INFO['ExAC_ALL'][0]) +"\t" +
# 		str(record.INFO['CLINSIG'][0]) +"\t" +
# 		str(record.INFO['Func.refGene'][0]) +"\t" +
# 		str(record.INFO['ExonicFunc.refGene'][0]) + "\t" +
# 		str(record.INFO['AAChange.refGene']) +"\t" +
# 		str(record.CHROM) + "\t" +
# 		str(record.POS) + "\t" +
# 		str(record.REF) + "\t" +
# 		str(record.ALT) + "\t" +
# 		str(record.INFO['dbscSNV_ADA_SCORE'][0]) + "\t" +
# 		str(record.INFO['dbscSNV_RF_SCORE'][0]) + "\t" +
# 		str(record.INFO['dpsi_zscore'][0]) + "\t" +
# 		str(s[0]) + "\t" +
# 		str(s[2]) + "\n"
# 	)
# f.close()
#
# os.system("sort -k1,1 -k2rn "+ outFile +  "> t.xls")
# os.system("echo \""+ head + "\" > " + outFile)
# os.system("cat t.xls >> " + outFile + " && rm t.xls")

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
        sys.exit('Multi references on vcf at "' + str(record.CHROM) + ':' + str(record.POS) + '". It seems that your vcf not followed the specifications.')

    if (len(record.ALT) > 1):
        sys.exit('Multi allelic variant on vcf at "' + str(record.CHROM) + ':' + str(record.POS) + '". See documentation and provide a well processed vcf (split variants).')

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
        "deleterioous":deleterious
    }

def is_clinvar_pathogenic(clinsig):
    """
    @summary: Define if clinvar annotation predict this variant as pathogenic
    @param clinsig: [str] The clinvar annotation provided by the vcf
    @return: [bool] True if is pathogenic and no Benign; False in other cases
    """
    # No clinsig available
    if clinsig == None:
        return False

    # Test if "Pathogenic" or "Benign" match on clinsig
    match_pathogenic = re.search("pathogenic", clinsig, re.IGNORECASE)
    match_benign = re.search("benign", clinsig, re.IGNORECASE)

    # Determine if clinvar as no doubt about pathogenicity
    if(match_pathogenic and not match_benign):
        return True
    else:
        return False


def is_splice_impact(splices_scores, is_indel, funcRefGene):
    """
    @summary: Predict splicing effect of the variant
    @param splices_scores: [dict] The dictionnary of splicing scores
    @param is_indel: [bool] Boolean to define if variants is indel or not
    @param funcRefGene: [str] Annotation provided by refGene about the biological function
    @return: [bool] True if is splicing impact; False in other cases
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
    if(ADA_splice or RF_splice or Zscore_splice or home_splice):
        return True
    else:
        return False

def is_stop_impact(exonicFuncRefGene, is_indel, funcRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @return: [bool] True if is stop impact; False in other cases
    """
    if(exonicFuncRefGene == None):
        return False


    match_stoploss = re.search("stoploss", exonicFuncRefGene, re.IGNORECASE)
    match_stopgain = re.search("stopgain", exonicFuncRefGene, re.IGNORECASE)

    if(match_stopgain or match_stoploss):
        return True
    else:
        return False

def is_frameshift_impact(exonicFuncRefGene, is_indel, funcRefGene):
    """
    @summary: Predict stop codon effect of the variant
    @param exonicFuncRefGene: [str] The exonic function predicted by RefGene
    @param is_indel: [bool] Boolean to define if variants is indel or not
    @param funcRefGene: [str] Annotation provided by refGene about the biological function
    @return: [bool] True if is frameshift impact; False in other cases
    """
    if(exonicFuncRefGene == None or funcRefGene == None):
        return False

    match_frameshift = re.search("frameshift", exonicFuncRefGene, re.IGNORECASE)
    match_exonic = re.search("exonic", funcRefGene, re.IGNORECASE)

    if(match_frameshift or (is_indel and match_exonic)):
        return True
    else:
        return False

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
                log.error(str(e))
                return

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

            splices_scores = {
                "ADA": record.INFO['dbscSNV_ADA_SCORE'][0],
                "RF": record.INFO['dbscSNV_RF_SCORE'][0],
                "Zscore":record.INFO['dpsi_zscore'][0],
            }

            # Calculate adjusted score for each variants
            adjusted_score = calculate_adjusted_score(impacts_scores)

            # Determine if variant is well annotated with clinvar as deleterious
            clinvar_pathogenicity = is_clinvar_pathogenic(record.INFO['CLINSIG'][0])

            # Determine the impact on splicing
            splice_impact = is_splice_impact(splices_scores, record.is_indel, record.INFO['Func.refGene'][0])

            # Determine the stop impact
            stop_impact = is_stop_impact(record.INFO['ExonicFunc.refGene'][0])

            # Determine the frameshift impact
            frameshift_impact = is_framshift_impact(record.INFO['ExonicFunc.refGene'][0],  record.is_indel, record.INFO['Func.refGene'][0])
            # print(record.INFO['ExonicFunc.refGene'])






########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Annotate VCF with Mobidic Prioritization Algorithm score (MPA).")
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)

    group_input = parser.add_argument_group('Inputs') # Inputs
    group_input.add_argument('-i', '--input', required=True, help="The vcf file to annotate (format: VCF). This vcf must be annotate with annovar.")

    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output', required=True, help="The csv file corresponding to the vcf file enter on input. (format : CSV)")
    args = parser.parse_args()


    # Process
    logging.basicConfig(format='%(asctime)s - %(name)s [%(levelname)s] %(message)s')
    log = logging.getLogger("MPA_score")
    log.setLevel(args.logging_level)
    log.info("Start MPA annotation")
    log.info("Command: " + " ".join(sys.argv))
    process(args, log)
    log.info("End MPA annotation")
