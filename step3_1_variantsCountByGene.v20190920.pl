#module load Perl/5.22.1-foss-2016a
#all human genes
use IO::Zlib;
use List::MoreUtils qw(first_index);
open IN_GENE, "xgen-exome-research-panel-gene-list9e255a1532796e2eaa53ff00001c1b3c.txt";
@genes;
while($line_gene=<IN_GENE>){
	chomp($line_gene);
	push (@genes,$line_gene);
}


my %gene_counts;
my %allele_counts;
open IN_SPL, "PCAWG_Sampleinfo.txt";
while($line_spl=<IN_SPL>){
	chomp($line_spl);
	foreach(@genes){
		$gene_counts{$line_spl}{$_}=0;
		$allele_counts{$line_spl}{$_}=0;
	}
}


$input=$ARGV[0];
$label=$ARGV[1];
my @maf= glob $input;

my $n=0;
foreach(@maf){
	chomp($_);
	$sub_maf=$_;
	$inmaf = new IO::Zlib;
	$inmaf->open($sub_maf, "rb");
	print $sub_maf."\n";
	#open INMAF, $sub_maf;
	while($line1=<$inmaf>){
		chomp($line1);
		@variants=split("\t",$line1);
		$len=$line1;
		#print join("\t",@variants)."\n";
		if(($n eq 0 )&& ($line1=~m/^Chrom.+/)){
			#@sample=@variants[10..$len-1];
			#print join("\t",@variants)."\n";
			($gene_index) = first_index{$_ eq "Gene.refGene"}@variants;
			($spl_index) = first_index{$_ eq "Sample"}@variants;
			($zyg_index) = first_index{$_ eq "zyg"}@variants;
			 #print $spl_index."\t".$gene_index."\n";
             #print $variants[$spl_index]."\t".$variants[$gene_index]."\n";
			open OUT_COUNT, ">PCAWG_".$label."_GeneLevelCounts.txt";
			open OUT_ALLELE, ">PCAWG_".$label."_GeneLevelAllelesCounts.txt";
			$n++;
		}else{
			if(exists($gene_counts{$variants[$spl_index]})){
				if(exists($gene_counts{$variants[$spl_index]}{$variants[$gene_index]})){
					$gene_counts{$variants[$spl_index]}{$variants[$gene_index]}=$gene_counts{$variants[$spl_index]}{$variants[$gene_index]}+1;
				}else{
					print $variants[$gene_index]." not included in the GeneList\n";
				}
			}else{
				print $variants[$spl_index]." not included in the SampleList\n";
			}

			if(exists($allele_counts{$variants[$spl_index]})){
				if(exists($allele_counts{$variants[$spl_index]}{$variants[$gene_index]})){
					if($variants[$zyg_index] eq "hom"){
						$allele_counts{$variants[$spl_index]}{$variants[$gene_index]}=$allele_counts{$variants[$spl_index]}{$variants[$gene_index]}+2;
					}else{
						$allele_counts{$variants[$spl_index]}{$variants[$gene_index]}=$allele_counts{$variants[$spl_index]}{$variants[$gene_index]}+1;
					}
				}
			}
		}
	}
}

@samples=keys %gene_counts;
$len=@samples;
print $len."\n";
print OUT_COUNT "\t".join("\t",@samples)."\n";
print OUT_ALLELE "\t".join("\t",@samples)."\n";

foreach(@genes){
	my $g=$_;
	print OUT_COUNT $g; 
	print OUT_ALLELE $g;
	foreach(@samples){
		my $s=$_;
		print OUT_COUNT "\t".$gene_counts{$s}{$g}; 
		print OUT_ALLELE "\t".$allele_counts{$s}{$g}; 
	}
	print OUT_COUNT "\n";
	print OUT_ALLELE "\n"
}






