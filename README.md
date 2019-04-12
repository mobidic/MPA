# Welcome to the MoBiDiC Prioritizing Algorithm (MPA)

- [Welcome to the MoBiDiC Prioritizing Algorithm (MPA) !](#welcome-to-the-mobidic-prioritizing-algorithm-mpa-)
	- [Overview](#overview)
		- [Citing MPA](#citing-mpa)
		- [Input](#input)
		- [Output](#output)
	- [Installation](#installation)
		- [Requirements](#requirements)
	- [Quick start](#quick-start)
		- [Quick guide for Annovar](#quick-guide-for-annovar)
			- [Install Annovar](#install-annovar)
			- [Download all databases](#download-all-databases)
			- [Annotate a VCF](#annotate-a-vcf)

--------------------------------------------------------------------------------

![MPA](https://raw.githubusercontent.com/mobidic/MPA/master/doc/img/logo-MPA.png)

![license](https://img.shields.io/github/license/mobidic/mpa.svg)
![release](https://img.shields.io/github/release/mobidic/mpa.svg)
![pypiV](https://img.shields.io/pypi/v/mobidic-mpa.svg)
![pythonV](https://img.shields.io/pypi/v/mobidic-mpa.svg)
![pypiS](https://img.shields.io/pypi/status/mobidic-mpa.svg)

## Overview

The MPA is a prioritizing algorithm for Next Generation Sequencing molecular
diagnosis. We propose an open source and free for academic user workflow.

Variant ranking is made with a unique score that take into account curated
database, biological assumptions, splicing predictions and the sum of various
predictors for missense alterations. Annotations are made for exonic and
splicing variants up to +300nt.

We show the pertinence of our clinical diagnosis approach with an updated
evaluation of in silico prediction tools using DYSF, DMD, LMNA, NEB and TTN
variants from the human expert-feeded Universal Mutation Database [1] with
courtesy regards of curators for pathogenic variants and from the ExAc database
[2] to define the dataset of neutral variants.

MPA needs an annotated vcf by ANNOVAR and give as output an annotated vcf with MPA score & ranks.

![MPA diagram](https://raw.githubusercontent.com/mobidic/MPA/master/doc/img/MPA_diagram2.png)

\*PTC: Premature Truncation Codon : nonsense or frameshift

\**: intronic positions between -20 and +5

### Citing MPA

> **Yauy et al.** MPA, a free, accessible and efficient pipeline for SNV annotation and prioritization for NGS routine molecular diagnosis. **The Journal of Molecular Diagnostics (2018)** https://doi.org/10.1016/j.jmoldx.2018.03.009

### Input

The MPA uses, as input, an annotated VCF file with Annovar [3] and the following
databases :

- Curated database: ClinVar [4]
- Biological assumption : refGene [5]
- Splicing predicition : SpliceAI [6], dbscSNV [7]
- Missense prediction : dbNSFP [8]

> Note : Short tutorial to annotate your VCF with Annovar (cf. [Quick guide for Annovar](#quick-guide-for-annovar)).

> **Update April 2019: spliceAI annotations now replace spidex. Waiting for spliceAI to be included in ANNOVAR, Files for this dataset in the proper format are available upon request (hg19 or hg38).**

> Multi-allelic variants in vcf should be splitted to biallelic variants with bcftools norm.

```bash
bcftools norm -m - file.vcf > file_breakmulti.vcf
```

### Output

#### In a VCF format

VCF is annotated with multiples items : MPA_impact (Clinvar_pathogenicity, splice_impact, stop and frameshift_impact, missense_impact and unknown_impact), MPA_ranking (1 to 8), MPA_final_score (from 0 to 10) and details for the scoring as MPA_available (from 0 to 10 missense tools which annotate), MPA_deleterious (number of missense tools that annotate pathogenic), MPA_ajusted (normalize missense score from 0 to 10).

#### Ranking : from 1 to 7 and score

- 1 - 10 with clinvar_pathogenicity : Pathogenic variants reported on ClinVar
- 2 - 10 with stop or frameshift_impact : Premature Truncation Codon : nonsense or frameshift
- 3,4,5 - 10 with splicing_impact (ADA, RF, spliceAI) : Affecting splice variants predictions ranked by algorithm performance robustness and strength
- 6 - with splicing_impact (indel) - Indel in splicing regions (as there is no splicing predictions for this case)
- 7 - with missense_impact (10 to 0) : Missense variants scores
- 8 - with unknown_impact : Exonic variants with not clearly annotated ORFs and splicing variants not predicted pathogenic

#### With a simple interface (Captain ACHAB)

MPA is a part of [MobiDL](https://github.com/mobidic/MobiDL) captainAchab workflow. MPA is the core of ranking in our useful and simple interface to easily interpret NGS variants at a glance named Captain ACHAB.
Find more informations at [Captain ACHAB](https://github.com/mobidic/Captain-ACHAB)

--------------------------------------------------------------------------------

## Installation

### Requirements

* Python 3

### pip

```bash
python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mobidic-mpa
```

## Quick start

To run the MPA script, use this command line :

```bash
mpa -i path/to/input.vcf -o path/to/output.vcf
```

### Quick guide for Annovar

This algorithm introduce here need some basics annotation. We introduce here a
quick guide to annotate your VCF files with Annovar.

#### Install Annovar

Follow instruction to download Annovar at :
> [http://www.openbioinformatics.org/annovar/annovar_download_form.php](http://www.openbioinformatics.org/annovar/annovar_download_form.php)

Unpack the package by using this command :

```bash
tar xvfz annovar.latest.tar.gz
```

#### Download all databases

In Annovar folder, download all database needed with annotate_variation.pl:

```bash
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp33a  humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbscsnv11 humandb/
```

Deprecated: For Spidex database, follow instruction here :

> [http://www.openbioinformatics.org/annovar/spidex_download_form.php](http://www.openbioinformatics.org/annovar/spidex_download_form.php)

> **Update April 2019: spliceAI annotations now replace spidex. Waiting for spliceAI to be included in ANNOVAR, Files for this dataset in the proper format are available upon request (hg19 or hg38).**

#### Annotate a VCF

The following command line annotate a VCF file :

```bash
perl path/to/table_annovar.pl path/to/example.vcf humandb/ -buildver hg19 -out path/to/output/name -remove -protocol refGene,refGene,clinvar_20180603,dbnsfp33a,spliceai_filtered,dbscsnv11 -operation g,g,f,f,f,f -nastring . -vcfinput -otherinfo -arg '-splicing 20','-hgvs',,,,
```

--------------------------------------------------------------------------------

**Montpellier Bioinformatique pour le Diagnostique Clinique (MoBiDiC)**

*CHU de Montpellier*

France

![MoBiDiC](https://github.com/mobidic/MPA/tree/master/doc/img/logo-mobidic.png)

[Visit our website](https://neuro-2.iurc.montp.inserm.fr/mobidic/)

--------------------------------------------------------------------------------

1. Béroud, C. et al. UMD (Universal Mutation Database): 2005 update. *Hum. Mutat.* **26**, 184–191 (2005).
2. Lek, M. et al. Analysis of protein-coding genetic variation in 60,706 humans. *Nature* **536**, 285–291 (2016).
3. Wang, K., Li, M. & Hakonarson, H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. *Nucleic Acids Res.* **38**, e164–e164 (2010).
4. Landrum, M. J. et al. ClinVar: public archive of interpretations of clinically relevant variants. *Nucleic Acids Res.* **44**, D862–D868 (2015).
5. O’Leary, N. A. et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Res.* **44**, D733–45 (2016).
6. Jaganathan et al. Predicting Splicing from Primary Sequence with Deep Learning. *Cell* **176**, 535-548 (2019).
7. Jian, X., Boerwinkle, E. & Liu, X. In silico prediction of splice-altering single nucleotide variants in the human genome. *Nucleic Acids Res.* **42**, 13534–13544 (2014).
8. Liu, X., Wu, C., Li, C. & Boerwinkle, E. dbNSFP v3.0: A One-Stop Database of Functional Predictions and Annotations for Human Nonsynonymous and Splice-Site SNVs. *Hum. Mutat.* **37**, 235–241 (2016).
