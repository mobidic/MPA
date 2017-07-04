##!/usr/local/bin/ python3.5
# coding: utf-8
# myoscore.py aim to calculate a unique score based on prediction tools in order
# to prioritize variants with a unique score

# ==============================================================================
# IMPORT LIBRARIES

import vcf		# read vcf => PyVCF :https://pyvcf.readthedocs.io/en/latest/
import sys		# system command
import csv		# read, write csv
import re		# regex
import argparse	# for options
import os		# for options
import pprint	# pretty print

PATH_SCORE_SERVER = 'score/score/'

# ==============================================================================

def parse(args):
	"""Parse arguments
    Keyword arguments:
	    args -- command line arguments
	Return :
		args.input	-- input file
		argd.output	-- output file
    """
	# Initialise parser
	parser = argparse.ArgumentParser(
		description='Calculating the meta score on annotated vcf')

	# Add group for mandatory arguments
	mandatoryArgs = parser.add_argument_group('mandatory arguments')
	## input file
	mandatoryArgs.add_argument('-i','--input', metavar='file.vcf',
		type=lambda arg: is_valid_path(parser, arg), help='vcf file annotated',
		required=True)
	## output file
	mandatoryArgs.add_argument('-o','--output', metavar='file.csv',
		type=argparse.FileType('w'),
		help='generated .csv file with variants and score', required=True)


	args = parser.parse_args()
	return(
		args.input,
		args.output,
	)
# ==============================================================================

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
		return(score_adjusted, "10db", available, "clinvar")
	# Test impact on stop codon
	if (exonicFunction =="stopgain" or exonicFunction =="stoploss"):
		return(score_adjusted, "10sfs", available, "stop")
	# Test impact on frameshift
	if (exonicFunction == "frameshift deletion" or
		exonicFunction == "frameshift insertion"):
		return(score_adjusted, "10sfs", available, "frameshift")
	# Test impact on splice function
	if ((ref or alt == "-") and re.match("splicing",function)):
		return(score_adjusted, "10sp",available, "splice")
	# Test impact on splice function
	if(str(ADAscore) != "."):
		if(float(ADAscore) >= 0.6):
			return(score_adjusted, "10sp",available, "splice")
	if(str(RFscore) != "."):
		if(float(RFscore) >= 0.6):
			return(score_adjusted, "10sp",available, "splice")
	if(str(Zscore) != "."):
		if(float(Zscore) < -2):
			return(score_adjusted, "10sp", available, "splice")
	# Test if variant maps to multiple location, return "u" (UNKNOWN) if it does
	if(exonicFunction == "unknown"):
		return(score_adjusted, "u",available, "u")
	#return 'na' if no tools score are available
	#return deleterious/available * 10, "%s/%s" % (deleterious, available)
	else:
		return score_adjusted, score_adjusted, available, "na"

# ==============================================================================
# ==============================================================================

def is_valid_path(parser, tested_path):
    '''
    Test if the tested_file is a valid file.
    A valid file is an existing file.
    :param parser: The argument parser
    :type parser: ArgumentParser
    :param tested_file: Tested file path
    :type tested_file: str
    :returns: tested_file
    :rtype: str
    '''
    if not os.path.exists(tested_path):
        parser.error("The path %s does not exist!" % tested_path)
    else:
        return tested_path  # return the path


# ==============================================================================

# Assigning command line arguments

(inFile, outFile) = parse(sys.argv)

try :
	inFile_h = open(inFile, 'r')
	vcf_reader = vcf.Reader(inFile_h)
	f = outFile
except UnicodeDecodeError :
	print ("\n")
	print("###*** The input file format is invalid. Please use an \
	an annotated VCF ***###")
	print ("\n")
	sys.exit()

# ==============================================================================

# Writing the header
f.write("Score\tGene.refGene\t")

# Writing sample name in the header
for elt in vcf_reader.samples:
	f.write(elt + "\t")
# Writing the header (second part)
f.write("ExAc\tClinSig\tFunc.refGene\t\
ExonicFunc.refGene\tAAChange.refGene\t#CHROM\tPOS\tREF\tALT\tADA\tRF\t\
Spidex\tScore\tNumber of tools\n")

# ==============================================================================

def test_vcf_annotation(record, key, error):

	try:
		record.INFO[key] = ["." if elt is None else elt for elt in record.INFO[key]]
	except KeyError:
		print(
			"The input VCF file provided lacks the annotation by " +
			str(key)
		)
		error = True
	return error

# ==============================================================================

#Create a list of the keys used in PyVCF and required for the scoring
vcf_key = [
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

# Iterating through the input VCF
for record in vcf_reader:

# =======																========
# Testing if the input vcf is correctly annotated

	error = False
	try:
		record.REF 	= ["." if record.REF is None else record.REF]
	except KeyError:
		print("The input VCF file provided lacks the annotation by REF.")
		error=True

	try:
		record.ALT 	= ["." if elt is None else elt for elt in record.ALT]
	except KeyError:
		print("The input VCF file provided lacks the annotation by ALT.")
		error=True

	for elt in vcf_key:
		error = test_vcf_annotation(record, elt, error)

	if error==True :
		sys.exit("\n" + "Please refer yourself to the documentation in order to\
 process VCF in a valid format" + "\n")

# =======                                                              =========

# call score function, instanciate it to "s"

	s = score(
		record.INFO['ExonicFunc.refGene'][0],
		record.REF,record.ALT[0],
		record.INFO['Func.refGene'][0],
		record.INFO['dbscSNV_ADA_SCORE'][0],
		record.INFO['dbscSNV_RF_SCORE'][0],
		record.INFO['dpsi_zscore'][0],
		record.INFO['SIFT_pred'][0],
		record.INFO['Polyphen2_HDIV_pred'][0],
		record.INFO['Polyphen2_HVAR_pred'][0],
		record.INFO['LRT_pred'][0],
		record.INFO['MutationTaster_pred'][0],
		record.INFO['FATHMM_pred'][0],
		record.INFO['PROVEAN_pred'][0],
		record.INFO['fathmm-MKL_coding_pred'][0],
		record.INFO['MetaSVM_pred'][0],
		record.INFO['MetaLR_pred'][0],
		record.INFO['CLINSIG']
	)

# ==============================================================================

# Writing output csv
	geno	= ""
	for sample in record.samples:
			geno += ("'" + sample['GT'] + "\t")
	f.write(
		str(s[1]) + "\t" +
		str(record.INFO['Gene.refGene'][0]) +"\t" +
		str(geno) +
		str(record.INFO['ExAC_ALL'][0]) +"\t" +
		str(record.INFO['CLINSIG'][0]) +"\t" +
		str(record.INFO['Func.refGene'][0]) +"\t" +
		str(record.INFO['ExonicFunc.refGene'][0]) + "\t" +
		str(record.INFO['AAChange.refGene']) +"\t" +
		str(record.CHROM) + "\t" +
		str(record.POS) + "\t" +
		str(record.REF) + "\t" +
		str(record.ALT) + "\t" +
		str(record.INFO['dbscSNV_ADA_SCORE'][0]) + "\t" +
		str(record.INFO['dbscSNV_RF_SCORE'][0]) + "\t" +
		str(record.INFO['dpsi_zscore'][0]) + "\t" +
		str(s[0]) + "\t" +
		str(s[2]) + "\n"
	)
f.close()

################################################################################
