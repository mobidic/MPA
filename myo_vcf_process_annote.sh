## Fusion de VCF
cd
cd /Users/henri/vcftools/src/perl/VCFPanel

x=1
for i in *.vcf ; do mv "$i" "${x}_$i" 
x=$((x+1))
done

x=1
for i in *.vcf ; do mv "$i" "${x}.vcf" 
x=$((x+1))
done

# Indexer et compresser 
# Ajuster bgzip et tabix -p vcf au nombre de fichier .vcf à traiter
for vcf in *.vcf ; do bgzip -f "$vcf"
done

for tab in *.gz ; do tabix -p vcf -f "$tab"
done


cd
cd vcftools/src/perl
perl vcf-merge VCFPanel/*.vcf.gz | bgzip -c > VCFPanel/MERGE.vcf.gz

cd
cd /Users/henri/vcftools/src/perl/VCFPanel



tabix MERGE.vcf.gz
# split multiallelic sites into biallelic records if both SNPs and indels should be merged separately into two records, specify both;

bcftools norm -m-both -o MERGE_split.vcf MERGE.vcf.gz

#-f, --fasta-ref FILE, reference sequence. Supplying this option will turn on left-alignment and normalization,
# Left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; recover multiallelics from multiple rows. 
# Left-alignment and normalization will only be applied if the --fasta-ref option is supplied

bcftools norm -f /Users/henri/Bioinfo/genome_ref/hg19.fa -o MERGE_split_normalized.vcf MERGE_split.vcf

# Copier le VCF merge dans le dossier annovar correspondant 
cp /Users/henri/vcftools/src/perl/VCFPanel/MERGE_split_normalized.vcf /Users/henri/Bioinfo/annovar/VCFbrut


### Annotation ANNOVAR

# Acceder au dossier
cd
cd /Users/henri/Bioinfo/annovar

# Ligne de commande annotation du vcf par annovar

perl /Users/henri/Bioinfo/annovar/table_annovar.pl \
	/Users/henri/Bioinfo/annovar/VCFbrut/MERGE_split_normalized.vcf \
	humandb/ \
	-buildver hg19 \
	-out /Users/henri/Bioinfo/annovar/VCFannote/MERGE_split_normalized_anno \
	-remove \
	-protocol refGene,knownGene,ccdsGene,cytoBand,clinvar_20170130,popfreq_max_20150413,kaviar_20150923,exac03,1000g2014oct_all,snp138,dbnsfp33a,spidex,dbscsnv11,mcap \
	-operation g,g,g,r,f,f,f,f,f,f,f,f,f,f \
	-nastring . \
	-vcfinput 	