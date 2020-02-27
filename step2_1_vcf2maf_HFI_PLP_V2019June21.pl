#module load Perl/5.22.1-foss-2016a
#vcf file name 
use IO::Zlib;

$vcf = $ARGV[0];
print $vcf."\n";
$spl = $ARGV[1];

mkdir gHFI;
$datatype="./gHFI/PCAWG";  #$ARGV[1];
$version="MetaSVM_ClinVar_v20190621";
#wirte file
$out_h=$datatype."_".$version."_".$spl."_HFI_PLP.txt.gz";
$OUT_H = new IO::Zlib;
$OUT_H->open($out_h, "wb9");

$out_t=$datatype."_".$version."_".$spl."_Neutral_PLP.txt.gz";
$OUT_T = new IO::Zlib;
$OUT_T->open($out_t, "wb9");

$out_u=$datatype."_".$version."_".$spl."_Unpredict_PLP.txt.gz";
$OUT_U = new IO::Zlib;
$OUT_U->open($out_u, "wb9");

$out_s=$datatype."_".$version."_".$spl."_Synonymous_PLP.txt.gz";
$OUT_S = new IO::Zlib;
$OUT_S->open($out_s, "wb9");

#clinVar database
$clinVar = new IO::Zlib;
$clinVar->open("variant_summary_2019-05.txt.gz", "rb");

%hash_clinvar;
while($line_clinvar=<$clinVar>){
        chomp($line_clinvar);
        @elem_clinvar=split("\t",$line_clinvar);
        #Chromosome|Start_Position|Reference_Allele|Tumor_Seq_Allele2|Sample
        #$elem[8]=~m/(TCGA-.+-.+)-.+-.+-.+-.+/;
        $mut_clinvar=$elem_clinvar[18]."|".$elem_clinvar[19]."|".$elem_clinvar[21]."|".$elem_clinvar[22];
        unless(exists($hash_clinvar{$mut_clinvar})){
                $hash_clinvar{$mut_clinvar}{"Guidelines"}=$elem_clinvar[26];
                $hash_clinvar{$mut_clinvar}{"TestedInGTR"}=$elem_clinvar[27];
        }
}

my $n=0;

$vcfFile = new IO::Zlib;
$vcfFile->open($vcf, "rb");
	#important iterms to extract
	@info_array=qw(Gene.refGene ExonicFunc.refGene AAChange.refGene CLNSIG CLNDN exome_AF genome_AF exome_AF_nfe genome_AF_nfe exome_AF_fin genome_AF_fin exome_non_cancer_AF_popmax genome_non_cancer_AF_popmax MetaSVM_rankscore MetaSVM_pred MetaLR_rankscore MetaLR_pred);
	while($line=<$vcfFile>){
		chomp($line);
		if($n==0 && $line=~m/#CHROM/){
			@elem_nam=split("\t",$line);
			$spl_num=@elem_nam-1;
			@spl_nam=@elem_nam[9..$spl_num];
			my $header="Chrom\tPos\tRef\tALT\tQual\tFilter\t".join("\t",@info_array)."\tSample\tGT\tzyg\tClinVarACMG\tClinVarGTR\n"; 
			print $OUT_H $header;
			print $OUT_T $header;
			print $OUT_U $header;
			print $OUT_S $header;
			$n++;
		}

		if($line!~m/^#.+/){
			@elem=split("\t",$line);
			@info=split(";",$elem[7]);
			#convert variant annotation info to array
			my %info_elem;
			foreach (@info) {
				@tmp=split("=",$_);
				$info_elem{$tmp[0]}= $tmp[1];
			}
			
			my @info_value=();
			foreach (@info_array) {
				push @info_value, $info_elem{$_};
			}

			$mut_HFI=$elem[0]."|".$elem[1]."|".$elem[3]."|".$elem[4];

			my $i=0;
			my $j=0;

			foreach(@elem[9..$spl_num]){
				my $genotype=$_;
				my %qc_hash=depth_info($elem[8],$genotype);
				#if(($qc_hash{"zyg"} ne ".") && ($qc_hash{"DP"} ne ".") && ($qc_hash{"DP"} ne 0) && ($qc_hash{"AD"} ne 0)){
				if($qc_hash{"zyg"} ne "."){
					#if(($qc_hash{"AD"}/$qc_hash{"DP"} > 0.3) && ($qc_hash{"AD"} > 5) && ($qc_hash{"GQ"} > 40)){
						if(exists($hash_clinvar{$mut_HFI})){
							$clinVarInfoACMG=$hash_clinvar{$mut_HFI}{"Guidelines"};
							$clinVarInfoGTR=$hash_clinvar{$mut_HFI}{"TestedInGTR"};
						}else{
							$clinVarInfoACMG=NA;
							$clinVarInfoGTR=NA;
						}

						my $output=join("\t",@elem[0,1,3,4,5,6])."\t".join("\t",@info_value)."\t".$spl_nam[$i]."\t".$qc_hash{"GT"}."\t".$qc_hash{"zyg"}."\t".$clinVarInfoACMG."\t".$clinVarInfoGTR."\n";

						if($info_elem{"ExonicFunc.refGene"} eq "nonsynonymous_SNV"){
							#HFI
							if(($info_elem{"CLNSIG"} eq "Pathogenic") || ($info_elem{"CLNSIG"} =~m/^Pathogenic.+/)|| ($info_elem{"CLNSIG"} eq "Likely_pathogenic")|| ($info_elem{"CLNSIG"}=~m/^Likely_pathogenic.+/)){
								print $OUT_H $output;
							}else{
								if(($info_elem{"CLNSIG"} eq "Benign") || ($info_elem{"CLNSIG"} eq "Likely_benign") || ($info_elem{"CLNSIG"} =~ m/^Benign.+/) || ($info_elem{"CLNSIG"} =~ m/^Likely_benign.+/)){
									print $OUT_T $output;
								}else{
									if($info_elem{"MetaSVM_pred"} eq "D"){
										print $OUT_H $output;
									}else{
										if($info_elem{"MetaSVM_pred"} eq "T"){
											print $OUT_T $output;
										}else{
											if($j eq 0){
												print $OUT_U $output;
												$j++;
											}
										}
									}
								}
							}
						}

						if($info_elem{"ExonicFunc.refGene"} eq "synonymous_SNV"){
										print $OUT_S $output;	
						}
					#}
				}
				$i++;
			}
		}
	}


sub depth_info
{
	my($key,$value)=@_;

		@key=split(":",$key);
		@value=split(":",$value);
				
		#construct hash
		%info;
		$i=0;
		foreach(@key){
			$info{$_}=$value[$i];
			$i++;
		}

		if($info{"GT"} eq "1|1"){
			$zyg="hom";
		}else{
			if($info{"GT"} eq "0|1"){
				$zyg="het";
			}else{
				$zyg=".";
			}
		}


		#if($info{"AD"}=~m/(\d+)\s(\d+)/){#this is not \t , could be white space
				#@ad=split("\t",$info{"AD"});
		#		print "TRUE:\t".$1."\t".$2."\n";
		#		$ref_dp=$1;#$ad[0];
		#		$alt_dp=$2;#$ad[1];
		#}else{
		#	if($info{"AD"}=~m/(\d+),(\d+)/){
				#@ad=split(",",$info{"AD"});
		#		$ref_dp=$1;#$ad[0];
		#		$alt_dp=$2;#$ad[1];
		#	}else{
		#		if($info{"AD"}=~m/(\d+)/ || $info{"AD"} eq "."){
		#			$ref_dp=$info{"RD"};
		#			$alt_dp=$info{"AD"};
		#		}else{
		#			print "ERROR: Unexpected AD:RD format\n";
		#			print $info{"AD"}."\n";
		#		}
		#	}
		#}
	#$dp=$depth_info{"DP"};
	#$gq=$depth_info{"GQ"};
	%depth_info;
	$depth_info{"GT"}=$info{"GT"};
	$depth_info{"zyg"}=$zyg;
	#$depth_info{"RD"}=$ref_dp;
	#$depth_info{"AD"}=$alt_dp;
	#$depth_info{"DP"}=$info{"DP"};
	#$depth_info{"GQ"}=$info{"GQ"};
	return(%depth_info);
}

 $OUT_H ->close;
 $OUT_T ->close;
 $OUT_U ->close;
 $OUT_S ->close;
 