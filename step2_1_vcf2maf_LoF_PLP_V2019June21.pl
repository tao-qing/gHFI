#module load Perl/5.22.1-foss-2016a
use IO::Zlib;
#vcf file name 
my $vcf = $ARGV[0];#glob("*_lifted_over.anno.hg19_multianno.vcf");
print $vcf."\n";
my $spl = $ARGV[1];

mkdir gLoF;

my $dataset="./gLoF/PCAWG"; 
my $version="gnomAD_ClinVar_v20190621";
#wirte file
my $out_h=$dataset."_".$version."_".$spl."_HCLoF_PLP.txt.gz";
$OUT_H = new IO::Zlib;
$OUT_H->open($out_h, "wb9");

my $out_b=$dataset."_".$version."_".$spl."_BLoF_PLP.txt.gz";
$OUT_B = new IO::Zlib;
$OUT_B->open($out_b, "wb9");

my $out_u=$dataset."_".$version."_".$spl."_NonLoF_PLP.txt.gz";
$OUT_U = new IO::Zlib;
$OUT_U->open($out_u, "wb9");


#read database
open IN_GNOMEAD, "all_gnomeAD_HC_LoF.bed";
my %lofDB;

while($line=<IN_GNOMEAD>){
	chomp($line);
	@elem=split("\t",$line);
	$mutlofdb=$elem[0]."|".$elem[1]."|".$elem[2]."|".$elem[3];
	unless(exists($lofDB{$mutlofdb})){
		$lofDB{$mutlofdb}="HC";
	}
}

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
        #print $mut."\n";
}


my $n=0;
$vcfFile = new IO::Zlib;
$vcfFile->open($vcf, "rb");
	#important iterms to extract
	@info_array=qw(Gene.refGene ExonicFunc.refGene AAChange.refGene CLNSIG CLNDN exome_AF genome_AF exome_AF_nfe genome_AF_nfe exome_AF_fin genome_AF_fin exome_non_cancer_AF_popmax genome_non_cancer_AF_popmax);
	
	@type=qw(stopgain stoploss splicing frameshift_deletion frameshift_insertion nonframeshift_deletion nonframeshift_insertion);

	while($line=<$vcfFile>){
		chomp($line);

		if($n==0 && $line=~m/#CHROM/){
			@elem_nam=split("\t",$line);
			$spl_num=@elem_nam-1;
			@spl_nam=@elem_nam[9..$spl_num];
			my $header="Chrom\tPos\tRef\tALT\tQual\tFilter\t".join("\t",@info_array)."\tSample\tGT\tzyg\tAD\tDP\tGQ\tgnomAD\tClinVarACMG\tClinVarGTR\n"; 
			
			print $OUT_H $header;
            print $OUT_B $header;
            print $OUT_U $header;
			$n++;
		}

		if($line!~m/^#/){
			my @elem=split("\t",$line);
			my @info=split(";",$elem[7]);

			my %info_elem;
			foreach (@info) {
				@tmp=split("=",$_);
				$info_elem{$tmp[0]}= $tmp[1];
			}
			
			my @info_value=();
			foreach (@info_array) {
				push @info_value, $info_elem{$_};
			}

			#$info_elem{"ExonicFunc.refGene"}=~s/\s/_/g;
			#print $info_elem{"ExonicFunc.refGene"}."\n";

			my $mutLoF=$elem[0]."|".$elem[1]."|".$elem[3]."|".$elem[4];

			if(exists($hash_clinvar{$mutLoF})){
				$clinVarInfoACMG=$hash_clinvar{$mutLoF}{"Guidelines"};
				$clinVarInfoGTR=$hash_clinvar{$mutLoF}{"TestedInGTR"};
			}else{
				$clinVarInfoACMG=NA;
				$clinVarInfoGTR=NA;
			}

			if(exists($lofDB{$mutLoF})){
				$gnomADInfo=$lofDB{$mutLoF};
			}else{
				$gnomADInfo=NA;
			}
			


			my $i=0;
			my $m=0;
	        foreach(@elem[9..$spl_num]){
	        		$genotype=$_;
	                my %qc_hash=depth_info($elem[8],$genotype);
					if($qc_hash{"zyg"} ne "."){

							my $output=join("\t",@elem[0,1,3,4,5,6])."\t".join("\t",@info_value)."\t".$spl_nam[$i]."\t".$qc_hash{"GT"}."\t".$qc_hash{"zyg"}."\t".$gnomADInfo."\t".$clinVarInfoACMG."\t".$clinVarInfoGTR."\n";

						if($info_elem{"ExonicFunc.refGene"}=~/.+frameshift.+/){
							#for indel use crateria: #Pindel-unique calls (at least 30x coverage and 20% VAF). according to kuan et al. Cell 2018
							#if(($qc_hash{"AD"}/$qc_hash{"DP"} > 0.2) && ($qc_hash{"DP"} > 30) && ($qc_hash{"GQ"} > 40)){

								if((grep /^$info_elem{"ExonicFunc.refGene"}$/, @type) || ($info_elem{"Func.refGene"} eq "splicing")){
									if((exists($lofDB{$mutLoF}) && ($info_elem{"CLNSIG"} ne "Benign")  && ($info_elem{"CLNSIG"} !~ m/^Benign.+/)  && ($info_elem{"CLNSIG"} ne "Likely_benign")  && ($info_elem{"CLNSIG"} !~ m/^Likely_benign.+/)) || ($info_elem{"CLNSIG"} eq "Pathogenic") || ($info_elem{"CLNSIG"} =~ m/^Pathogenic.+/) || ($info_elem{"CLNSIG"} eq "Likely_pathogenic")|| ($info_elem{"CLNSIG"}=~m/^Likely_pathogenic.+/)){#high condidence or pathogenic LoF
										print $OUT_H $output;
									}else{
										if(($info_elem{"CLNSIG"} eq "Benign")||($info_elem{"CLNSIG"} =~ m/^Benign.+/)||($info_elem{"CLNSIG"} eq "Likely_benign")||($info_elem{"CLNSIG"} =~ m/^Likely_benign.+/)){
											print $OUT_B $output;
										}else{
											if($m eq 0){
												print $OUT_U $output;
												$m++;
											}
										}
									}
								}
			                #}
			            }else{
			            	#if(($qc_hash{"AD"}/$qc_hash{"DP"} > 0.3) && ($qc_hash{"AD"} > 5) && ($qc_hash{"GQ"} > 40)){
								if((grep /^$info_elem{"ExonicFunc.refGene"}$/, @type) || ($info_elem{"Func.refGene"} eq "splicing")){

									if((exists($lofDB{$mutLoF}) && ($info_elem{"CLNSIG"} ne "Benign")  && ($info_elem{"CLNSIG"} !~ m/^Benign.+/)  && ($info_elem{"CLNSIG"} ne "Likely_benign")  && ($info_elem{"CLNSIG"} !~ m/^Likely_benign.+/)) || ($info_elem{"CLNSIG"} eq "Pathogenic") || ($info_elem{"CLNSIG"} =~ m/^Pathogenic.+/)){#high condidence or pathogenic LoF
										print $OUT_H $output;
									}else{
										if(($info_elem{"CLNSIG"} eq "Benign")||($info_elem{"CLNSIG"} =~ m/^Benign.+/)||($info_elem{"CLNSIG"} eq "Likely_benign")||($info_elem{"CLNSIG"} =~ m/^Likely_benign.+/)){
											print $OUT_B $output;
										}else{
											if($m eq 0){
												print $OUT_U $output;
												$m++;
											}
										}
									}

								}
			                #}
			            }
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
		my %info;
		my $i=0;
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
		#		#@ad=split(",",$info{"AD"});
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
 $OUT_B ->close;
 $OUT_U ->close;

 