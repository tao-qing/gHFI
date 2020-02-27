#this workflow is to calculate the germline high-functional impact variants load for individuals.


####the dependent files####
ClinVar Database: variant_summary_2019-05.txt.gz 
gnomAD HCLoF: all_gnomeAD_HC_LoF.bed
All human genes: xgen-exome-research-panel-gene-list9e255a1532796e2eaa53ff00001c1b3c.txt


sample names, each column is one sample, the sample name has to be the same as the sample name in the vcf file.
Sampleinfo.txt

Put all these files under working directory. 

####1. Annotation####
annotate vcf file using Annovar.
need the follow annovar database:
dbnsfp30a, exac03, gnomad_exome, avsnp150, clinvar_20190305, cosmic70

those database can be installed: 
perl annotate_variation.pl -webfrom annovar -downdb [*database name*] -buildver hg19 

when all the annovar databases are ready.
run the annotation script: step1_1_AnnovarAnnotation.sh
change those file path accordingly.


####2. HFI/LoF identification####
define two kinds of HFI variants: Deleterious Missense, Loss-of-function.

Run the follow two script to get the HFI variants

module load Perl/5.22.1-foss-2016a (this script depends on some perl modules, use this version of perl)
perl step2_1_vcf2maf_HFI_PLP_V2019June21.pl [annotated vcf files from Step 1]  gHFI
perl step2_1_vcf2maf_LoF_PLP_V2019June21.pl [annotated vcf files from Step 1]  gLoF

####3. Summarize the HFI count####
In this step, it will calculate how many HFI/LoF variants for each huaman genes of each individuals. 

module load Perl/5.22.1-foss-2016a 
perl step3_1_variantsCountByGene.v20190920.pl *_HFI_PLP.txt.gz gHFI

perl step3_1_variantsCountByGene.v20190920.pl *_HCLoF_PLP.txt.gz gLoF
