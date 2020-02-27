#!/bin/bash
#SBATCH --array=0
#SBATCH --partition=bigmem,scavenge
#SBATCH --job-name=Plink
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --mem=200012MB

/home/tq37/software/annovar/table_annovar.pl -vcfinput *.vcf.gz /home/tq37/software/annovar/humandb/ -buildver hg19 -out *.anno -remove -protocol refGene,avsnp150,dbnsfp33a,gnomad211_exome,gnomad211_genome,clinvar_20190305 -operation gx,f,f,f,f,f -nastring . --polish; gzip -c  *.anno.hg19_multianno.vcf > *.anno.hg19_multianno.vcf.gz