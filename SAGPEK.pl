#!/usr/bin/perl -w
use strict;
use IO::File;
use lib "./lib";
use ABI;
use GD;
use Getopt::Long;
my ($type,$orientation,$tag_file,$AA,$chromatogram);


die "Wrong: must provide at least one parameter!\n" unless @ARGV;


GetOptions(
	"type=s" => \$type,
	"orientation:s" => \$orientation,
	"tag:s" => \$tag_file,
	"AA:s" => \$AA,
	"chromatogram:s" => \$chromatogram,

) or die "Usage: $0 [-type=type of analysis] [-orientation=r or f] [-tag=Tag file of Custom type] [-AA=on or off] [-chromatogram=on or off] \n";

$orientation //= "f";
$AA ="on" if defined $AA && $AA eq '';
$chromatogram = "off" if defined $chromatogram && $chromatogram eq '';
print "\nParameter Information:\n\tThe Orientation of Sanger sequencing is: $orientation!\n\tThe 'r' is for Reverse and 'f' for Forward!\n" if  $orientation;
#print "\tThe AA report function is: $AA!\n";
print "------------------------------------------------------------------------------\n";

#my $type = shift @ARGV;
#my $type="TEST";
#my $orientation = "f";

my %custom;

if($type=~/custom/i){
	if(defined $tag_file){
		print ("\tThe tag file for the Custom type is: $tag_file!\n");
		unless (-e "Custom/$tag_file") {
		die("Error Information:\n\tTag file $tag_file not found in the Custom directory!\n");
		}
	}else{
		die ("Error Information:\n\tTag file must be defined if you use the Custom type!\n");
	}
	if($AA=~/off|on/i){
		print ("\tThe Amino Acids report function is: $AA!\n");
	}else{
		die ("Error Information:\n\tThe AA function must be: on or off!\n");
	}
	print "------------------------------------------------------------------------------\n";
	my $tag_fh = IO::File->new("Custom/$tag_file",'r');
	#print "type is custom\ntag file is $tag_file\n";
	my $custom_tag_count=0;
	while(<$tag_fh>){
		chomp;
		my $line = $_;
		$custom_tag_count++;
		my ($tag_name, $tag_seq)=split /\s+/, $line;
		$custom{$custom_tag_count}{tag_name}=$tag_name;
		$custom{$custom_tag_count}{tag_seq}=$tag_seq;
	}
}

my @abifiles=<ABI/*ab1>;

if(@abifiles <1){
	die(" Error: No ABI files found! Please make sure ABI-format files are listed under the ABI directory!\n");
}

my $fh_geno_out;

if($type=~/^TEST$/i){
	$fh_geno_out=IO::File->new(">Genotype/TEST.genotype.txt");
	$fh_geno_out->print("SampleID");
	if($AA=~/on/i){
		$fh_geno_out->print("\tFirst_Site_Genotype\tFirst_site_AA\tSecond_Site_Genotype\tSecond_site_AA\tThird_Site_Genotype\tThird_Site_AA\n");
	}elsif($AA=~/off/i){
		$fh_geno_out->print("\tFirst_Site_Genotype\tSecond_Site_Genotype\tThird_Site_Genotype\n");
	}
}elsif($type=~/^cattle_CN$/i){
	$fh_geno_out= IO::File->new(">Genotype/cattle_CN.genotype.txt");
	$fh_geno_out->print("SampleID\tCN_genotype\tAmino Acids\n");
}elsif($type=~/^cattle_DUMPS$/){
	$fh_geno_out=IO::File->new(">Genotype/cattle_DUMPS.genotype.txt");
	$fh_geno_out->print("SampleID\tDUMPS_genotype\tAmino Acids\n");
}elsif($type=~/^HBV$/i){
	$fh_geno_out=IO::File->new(">Genotype/HBV.genotype.txt");
	$fh_geno_out->print("SampleID");
	my @rt_nums=(166, 169, 173, 180, 181, 184, 191, 194, 200, 202, 204, 207, 213, 214, 215, 217, 221, 229, 233, 236, 237, 238, 245, 250, 256);
	for my $rt_num(@rt_nums) {
		$fh_geno_out->print("\trt$rt_num");
	}
	$fh_geno_out->print("\n");
}elsif($type=~/^PAH$/){
	$fh_geno_out = IO::File->new(">Genotype/Newborn_PAH.genotype.txt");
	$fh_geno_out->print("SampleID");
	my @AA_nums = (53, 107, 111, 171, 223, 225, 243, 304, 356, 388, 392, 399, 408, 413);
	for my $AA_num (@AA_nums) {
		$fh_geno_out->print("\tAA_${AA_num}_condons\tAA_${AA_num}_Amino_Acids");
	}
	$fh_geno_out->print("\n");
}elsif($type=~/^custom$/i){
	$fh_geno_out = IO::File->new(">Genotype/Custom.genotype.txt");
	$fh_geno_out ->print("SampleID");
	for my $custom_tag_count (sort {$a<=>$b} keys %custom) {
		my $custom_tag_name = $custom{$custom_tag_count}{tag_name};
		if($AA=~/on/i){
			$fh_geno_out->print("\t$custom_tag_name-condons\t$custom_tag_name-AA");
		}elsif($AA=~/off/i){
			$fh_geno_out->print("\t$custom_tag_name-condons");
		}
	}
	$fh_geno_out->print("\n");
}else{
	die "Error Information: the type \"$type\" is not supported by SAGPEK!\n";
}

for my $abi_file (@abifiles) {
	my $prefix;
	if($abi_file=~/ABI\/(.+)\.ab1/){
		$prefix=$1;
	}
	print("Processing: $prefix!\n");
	print("  Step 1: Extracting signals from ABI-format files......");
	
	my $amppeaks_result=getAMPpeaks($abi_file);
	print("Done!\n");

	print("  Step 2: Performing genotyping......");

	$fh_geno_out->print("$prefix");
	my $primary_seq;
	my $secondary_seq;
	my @sigs;
	my %complementary = (
		"A" => "T",
		"T" => "A",
		"C" => "G",
		"G" => "C"
	);
	
	my @amppeaks=split /\n/, $amppeaks_result;
	for my $ampsignal(@amppeaks) {
		my %data;
		my ($A,$C,$G,$T,$ratio,$sig)=split /\t/, $ampsignal;
		my %bases_index=(0=>"A",1=>"C",2=>"G",3=>"T");
		my @unsorted=($A,$C,$G,$T);
		my @sorted_indexes=sort {$unsorted[$b]<=>$unsorted[$a]} 0..$#unsorted;
		my @sorted_values=@unsorted[@sorted_indexes];
		my ($primary_index,$secondary_index)=@sorted_indexes[0,1];
		my $primary_base=$bases_index{$primary_index};
		$primary_seq.=$primary_base;
		my $secondary_base=$bases_index{$secondary_index};
		$secondary_seq.=$secondary_base;
		push @sigs, $sig;
	}

	if($orientation=~/^r$/i){
		my @primary= split //, $primary_seq;
		my @primary_r = reverse @primary;
		my @secondary = split //, $secondary_seq;
		my @secondary_r = reverse @secondary;
		my $primary_r_c="";
		my $secondary_r_c="";
		for my $p_r (@primary_r){
			my $p_r_c = $complementary{$p_r};
			$primary_r_c.=$p_r_c;
		}
		for my $s_r (@secondary_r){
			my $s_r_c = $complementary{$s_r};
			$secondary_r_c.=$s_r_c;
		}

		$primary_seq=$primary_r_c;
		$secondary_seq = $secondary_r_c;
	}

	my %condons = (
		"TTT" => "F",
		"TTC" => "F",
		"TTA" => "L",
		"TTG" => "L",
		"CTT" => "L",
		"CTC" => "L",
		"CTA" => "L",
		"CTG" => "L",
		"ATT" => "I",
		"ATC" => "I",
		"ATA" => "I",
		"ATG" => "M",
		"GTT" => "V",
		"GTC" => "V",
		"GTA" => "V",
		"GTG" => "V",
		"TCT" => "S",
		"TCC" => "S",
		"TCA" => "S",
		"TCG" => "S",
		"CCT" => "P",
		"CCC" => "P",
		"CCA" => "P",
		"CCG" => "P",
		"ACT" => "T",
		"ACC" => "T",
		"ACA" => "T",
		"ACG" => "T",
		"GCT" => "A",
		"GCC" => "A",
		"GCA" => "A",
		"GCG" => "A",
		"TAT" => "Y",
		"TAC" => "Y",
		"CAT" => "H",
		"CAC" => "H",
		"CAA" => "Q",
		"CAG" => "Q",
		"AAT" => "N",
		"AAC" => "N",
		"AAA" => "K",
		"AAG" => "K",
		"GAT" => "D",
		"GAC" => "D",
		"GAA" => "E",
		"GAG" => "E",
		"TGT" => "C",
		"TGC" => "C",
		"TGG" => "W",
		"CGT" => "R",
		"CGC" => "R",
		"CGA" => "R",
		"CGG" => "R",
		"AGT" => "S",
		"AGC" => "S",
		"AGA" => "R",
		"AGG" => "R",
		"GGT" => "G",
		"GGC" => "G",
		"GGA" => "G",
		"GGG" => "G",
		"TAA" => "STOP",
		"TAG" => "STOP",
		"TGA" => "STOP"
	);

	#####...#####
	# For TEST.
	if($type=~/^TEST$/i){	
		my %sites_based_on_1 = map_TEST_sites($primary_seq);
		for my $site_no (1, 2, 3) {
			my $tag_name="TEST_$site_no";
			my $condon_all_possible_string="NA";
			my $AA_all_possible_string="NA";
			if(exists $sites_based_on_1{$tag_name}){
				my $site_based_on_1 = $sites_based_on_1{$tag_name};

				my %condon_base;
				for my $condon_site (1..3) {
					my $site_index = $site_based_on_1- 2 + $condon_site;
					my $primary_base = substr($primary_seq, $site_index, 1);
					my $secondary_base = substr($secondary_seq, $site_index, 1);
					my $sig = $sigs[$site_index];
					if($sig eq "TRUE") {
						$condon_base{"base$condon_site"}{$primary_base}=1;
						$condon_base{"base$condon_site"}{$secondary_base}=1;
					}else{
						$condon_base{"base$condon_site"}{$primary_base}=1;
					}
				}
				my @condon_all_possible;
				for my $condon_1_base (keys %{$condon_base{base1}}) {
					for my $condon_2_base (keys %{$condon_base{base2}} ) {
						for my $condon_3_base (keys %{$condon_base{base3}}) {
							my $condon_zuhe="$condon_1_base$condon_2_base$condon_3_base";
							push @condon_all_possible, $condon_zuhe;
						}
					}

				}
				$condon_all_possible_string=join(",", @condon_all_possible);

				my %AA_all_possible;
				for my $condon_anyone (@condon_all_possible) {
					my $AA_anyone = $condons{$condon_anyone};
					$AA_all_possible{$AA_anyone}=1;
				}
				my @AA_all_possible=keys %AA_all_possible;
				$AA_all_possible_string=join(",", @AA_all_possible);
			}
			if($AA=~/^on$/i){
				$fh_geno_out->print("\t$condon_all_possible_string\t$AA_all_possible_string");
			}elsif($AA=~/^off$/i){
				$fh_geno_out->print("\t$condon_all_possible_string");
			}
		}
		$fh_geno_out->print("\n");
	}


	#####...#####
	# For HBV.
	if($type=~/^HBV$/i){
		my %HBV_ref=(
			"rt166"=>"F",
			"rt169"=>"I",
			"rt173"=>"V",
			"rt180"=>"L",
			"rt181"=>"A",
			"rt184"=>"T",
			"rt191"=>"V",
			"rt194"=>"A",
			"rt200"=>"A",
			"rt202"=>"S",
			"rt204"=>"M",
			"rt207"=>"V",
			"rt213"=>"S",
			"rt214"=>"V",
			"rt215"=>"Q",
			"rt217"=>"R",
			"rt221"=>"Y",
			"rt229"=>"L",
			"rt233"=>"I",
			"rt236"=>"N",
			"rt237"=>"P",
			"rt238"=>"N",
			"rt245"=>"Y",
			"rt250"=>"M",
			"rt256"=>"S"
		);

		my %sites_cor_based_on_1 = map_HBV_sites($primary_seq);
		my @rt_nums=(166, 169, 173, 180, 181, 184, 191, 194, 200, 202, 204, 207, 213, 214, 215, 217, 221, 229, 233, 236, 237, 238, 245, 250, 256);
		for my $rt_AA_num(@rt_nums) {
			if(exists $sites_cor_based_on_1{"rt$rt_AA_num"}){
				
				my $site_based_on_1 = $sites_cor_based_on_1{"rt$rt_AA_num"};
				my $site_index=$site_based_on_1-1;
				my $rt_condon=substr($primary_seq, $site_index, 3);
				my $rt_ref=$HBV_ref{"rt$rt_AA_num"};
				my $rt_detect=$condons{$rt_condon};
				if($rt_ref eq $rt_detect){
					$fh_geno_out->print("\t$rt_ref->$rt_detect.No mutant!");
				}else{
					$fh_geno_out->print("\t$rt_ref->$rt_detect");
				}
			}
			else{
				$fh_geno_out->print("\tNA");
			}
		}
		$fh_geno_out->print("\n");
	}

	#####...#####
	# For cattle.
	if($type =~/cattle/i){

		my %cattle_tags=(
			"cattle_CN"=> "TGTACGAGGAC",
			"cattle_DUMPS" => "TTCTGGCTCC"
			);
		my $cattle_tag=$cattle_tags{$type};
		my @sites = map_cattle_sites($primary_seq, $cattle_tag);
		my $geno="NA";
		my $AA="NA";
		if(@sites>0){
			my $site_based_on_1 = shift @sites;
			my $site_index = $site_based_on_1-1;
			my $primary_BASE = substr($primary_seq, $site_index, 1);
			my $secondary_BASE = substr($secondary_seq, $site_index, 1);
			my $sig=$sigs[$site_index];
			my $primary_condon="${primary_BASE}GA";
			my $secondary_condon="${secondary_BASE}GA";
			my $primary_AA=$condons{$primary_condon};
			my $secondary_AA=$condons{$secondary_condon};
			if($sig eq "TRUE"){
				$fh_geno_out->print("\t$primary_BASE/$secondary_BASE\t$primary_AA/$secondary_AA\n");
			}
			else{
				$fh_geno_out->print("\t$primary_BASE/$primary_BASE\t$primary_AA/$primary_AA\n");
			}
		}
		else{
			$fh_geno_out->print("\t$geno\t$AA\n");
		}
	}

	###...###
	# For newborn.
	if($type =~/^PAH$/i){
		my %sites_based_on_1 = map_newborn_PAH_sites ($primary_seq);
		my @AA_sites=(48, 65, 241, 243, 245, 261, 300, 388, 390, 403, 408, 413, 414, 415);
		for my $AA_site (@AA_sites) {
			my $condon_all_possible_string="NA";
			my $AA_all_possible_string="NA";
			if(exists $sites_based_on_1{"PAH_$AA_site"}){

				my $site_based_on_1 = $sites_based_on_1{"PAH_$AA_site"};

				my %target_condon;
				for my $condon_site (1..3) {
					my $site_index = $site_based_on_1 -2 + $condon_site;
					my $primary_base = substr($primary_seq, $site_index, 1);
					my $secondary_base = substr($secondary_seq, $site_index, 1);
					my $sig = $sigs[$site_index];
					if($sig eq "TRUE") {
						$target_condon{"condon_site_$condon_site"}{$primary_base}=1;
						$target_condon{"condon_site_$condon_site"}{$secondary_base}=1;
					}else{
						$target_condon{"condon_site_$condon_site"}{$primary_base}=1;
					}
				}
				my @condon_all_possible;
				for my $condon_1_base (keys %{$target_condon{condon_site_1}}) {
					for my $condon_2_base (keys %{$target_condon{condon_site_2}} ) {
						for my $condon_3_base (keys %{$target_condon{condon_site_3}}) {
							my $condon_zuhe="$condon_1_base$condon_2_base$condon_3_base";
							push @condon_all_possible, $condon_zuhe;
						}
					}

				}
				$condon_all_possible_string=join(",", @condon_all_possible);

				my %AA_all_possible;
				for my $condon_anyone (@condon_all_possible) {
					my $AA_anyone = $condons{$condon_anyone};
					$AA_all_possible{$AA_anyone}=1;
				}
				my @AA_all_possible=keys %AA_all_possible;
				$AA_all_possible_string=join(",", @AA_all_possible);
			}
			$fh_geno_out->print("\t$condon_all_possible_string\t$AA_all_possible_string");
		}
		$fh_geno_out->print("\n");
	}


	###...###
	# For Custom Sites.
	if($type =~/^custom$/i){
		my %sites_based_on_1 = map_Custom_sites ($primary_seq, \%custom);
		for my $custom_tag_count (sort {$a<=>$b} keys %custom) {
			my $tag_name = $custom{$custom_tag_count}{tag_name};
			my $tag_seq = $custom{$custom_tag_count}{tag_seq};
			my $condon_all_possible_string="NA";
			my $AA_all_possible_string="NA";
			if(exists $sites_based_on_1{$tag_name}){
				my $site_based_on_1 = $sites_based_on_1{$tag_name};

				my %target_condon;
				for my $condon_site (1..3) {
					my $site_index = $site_based_on_1 -2 + $condon_site;
					my $primary_base = substr($primary_seq, $site_index, 1);
					my $secondary_base = substr($secondary_seq, $site_index, 1);
					my $sig = $sigs[$site_index];
					if($sig eq "TRUE") {
						$target_condon{"condon_site_$condon_site"}{$primary_base}=1;
						$target_condon{"condon_site_$condon_site"}{$secondary_base}=1;
					}else{
						$target_condon{"condon_site_$condon_site"}{$primary_base}=1;
					}
				}
				my @condon_all_possible;
				for my $condon_1_base (keys %{$target_condon{condon_site_1}}) {
					for my $condon_2_base (keys %{$target_condon{condon_site_2}} ) {
						for my $condon_3_base (keys %{$target_condon{condon_site_3}}) {
							my $condon_zuhe="$condon_1_base$condon_2_base$condon_3_base";
							push @condon_all_possible, $condon_zuhe;
						}
					}

				}
				$condon_all_possible_string=join(",", @condon_all_possible);
				if(@condon_all_possible>2){
					$condon_all_possible_string="Warning: Two situations happen here: 1) two variations co-exit in the codon. 2) To much noise exist for this Sanger sequencing result!";
				}

				my %AA_all_possible;
				for my $condon_anyone (@condon_all_possible) {
					my $AA_anyone = $condons{$condon_anyone};
					$AA_all_possible{$AA_anyone}=1;
				}
				my @AA_all_possible=keys %AA_all_possible;
				
				$AA_all_possible_string=join(",", @AA_all_possible);
				if(@AA_all_possible>2){
					$AA_all_possible_string=".";
				}elsif(@AA_all_possible <1){
					$AA_all_possible_string="NA";
				}
			}
			if($AA=~/on/i){
				$fh_geno_out->print("\t$condon_all_possible_string\t$AA_all_possible_string");
			}elsif($AA=~/off/i){
				$fh_geno_out->print("\t$condon_all_possible_string");
			}
		}
		$fh_geno_out->print("\n");
	}
	if ($chromatogram=~/^on$/i) {
		plotchroma($abi_file);
	}

	print "Done!\n";
}

sub map_cattle_sites {
	my $seq = shift;
	my $cattle_tag = shift;
	my @sites_cor_based_on_1;
	if($seq=~/$cattle_tag/){
		my $previous=$`;
		my $matched=$&;
		my $post=$';
		my $seq_before_site="$previous$matched";
		my $coor_site=length($seq_before_site)+1;
		#print "$tag: $coor_site\n";
		push @sites_cor_based_on_1, $coor_site;
	}
	else{
		my @tag_bases=split //, $cattle_tag;
		my $match_flag=0;
		for my $base_pos (0..$#tag_bases) {
			my @new_tag_bases=@tag_bases;
			$new_tag_bases[$base_pos]="[ATGC]";
			
			my $new_tag=join("", @new_tag_bases);
			if($seq=~/$new_tag/){
				my $previous=$`;
				my $matched=$&;
				my $post=$';
				my $seq_before_site="$previous$matched";
				my $coor_site=length($seq_before_site)+1;
				#print "$new_tag: $coor_site\n";
				push @sites_cor_based_on_1, $coor_site;
				$match_flag++;
			}
		}
	}
	return @sites_cor_based_on_1;
}



sub map_TEST_sites {
	my $seq=shift;
	my %sites_cor_based_on_1;

	my %tags=(
		"TEST_1" => "TGGGCCCATC",
		"TEST_2" => "GCCGCCTTTC",
		"TEST_3" => "GCCTGAAGTA"
	);


	for my $tag_name (keys %tags) {
		my $tag = $tags{$tag_name};
		if($seq=~/$tag/){
			my $previous=$`;
			my $matched=$&;
			my $post=$';
			my $seq_before_site="$previous$matched";
			my $coor_site=length($seq_before_site)+1;
			$sites_cor_based_on_1{$tag_name}=$coor_site;
		}
		else{
			my @tag_bases=split //, $tag;
			my $match_flag=0;
			for my $base_pos (0..$#tag_bases) {
				my @new_tag_bases=@tag_bases;
				$new_tag_bases[$base_pos]="[ATGC]";
				
				my $new_tag=join("", @new_tag_bases);
				if($seq=~/$new_tag/){
					my $previous=$`;
					my $matched=$&;
					my $post=$';
					my $seq_before_site="$previous$matched";
					my $coor_site=length($seq_before_site)+1;
					$sites_cor_based_on_1{$tag_name}=$coor_site;
				}
			}
		}
	}
	return %sites_cor_based_on_1;
}


sub map_newborn_PAH_sites {
	my $seq = shift;
	my %tags=(
		"PAH_48" => "AGTTGGTGCA",
		"PAH_65" => "CCTGACCCAC",
		"PAH_241" => "CACTGGTTTC",
		"PAH_243" => "TTTCCGCCTC",
		"PAH_245" => "CCGCCGACCT",
		"PAH_261" => "CCTGGCCTTC",
		"PAH_300" => "TCGCAGCTTT",
		"PAH_388" => "CCTCTATTAC",
		"PAH_390" => "TTACGTGGCA",
		"PAH_403" => "AAGGAACTTT",
		"PAH_408" => "CACAATACCT",
		"PAH_413" => "CTTCTCAGTT", 
		"PAH_414" => "CTCAGTTCGC",
		"PAH_415" => "AGTTCGCTAC"
 	);
	my %sites_cor_based_on_1;
	for my $tag_name (keys %tags) {
		my $tag=$tags{$tag_name};

		if($tag_name=~/PAH_(\d+)/){
			my $tag_num=$1;
			my $tag=$tags{$tag_name};
			if($seq=~/$tag/){
				my $previous=$`;
				my $matched=$&;
				my $post=$';
				my $seq_before_site="$previous$matched";
				my $coor_site=length($seq_before_site)+1;
				$sites_cor_based_on_1{"PAH_$tag_num"}=$coor_site;
			}
			else{
				my @tag_bases=split //, $tag;
				my $match_flag=0;
				for my $base_pos (0..$#tag_bases) {
					my @new_tag_bases=@tag_bases;
					$new_tag_bases[$base_pos]="[ATGC]";
					
					my $new_tag=join("", @new_tag_bases);
					if($seq=~/$new_tag/){
						my $previous=$`;
						my $matched=$&;
						my $post=$';
						my $seq_before_site="$previous$matched";
						my $coor_site=length($seq_before_site)+1;
						$sites_cor_based_on_1{"PAH_$tag_num"}=$coor_site;
					}
				}
			}
		}
	}
	return %sites_cor_based_on_1;
}


sub map_Custom_sites {
	my $seq = shift;
	my $custom_ref = shift;
	my %custom1 = %$custom_ref;
	my %sites_cor_based_on_1;
	for my $custom_tag_count (sort {$a<=>$b} keys %custom1) {
		my $tag_name = $custom1{$custom_tag_count}{tag_name};
		my $tag_seq = $custom1{$custom_tag_count}{tag_seq};
		if($seq=~/$tag_seq/){
			my $previous=$`;
			my $matched=$&;
			my $post=$';
			my $seq_before_site="$previous$matched";
			my $coor_site=length($seq_before_site)+1;
			$sites_cor_based_on_1{$tag_name}=$coor_site;
		}
		else{
			my @tag_bases=split //, $tag_seq;
			my $match_flag=0;
			for my $base_pos (0..$#tag_bases) {
				my @new_tag_bases=@tag_bases;
				$new_tag_bases[$base_pos]="[ATGC]";
				my $new_tag=join("", @new_tag_bases);
				if($seq=~/$new_tag/){
					my $previous=$`;
					my $matched=$&;
					my $post=$';
					my $seq_before_site="$previous$matched";
					my $coor_site=length($seq_before_site)+1;
					$sites_cor_based_on_1{$tag_name}=$coor_site;
				}
			}
		}
	}
	return %sites_cor_based_on_1;
}


sub map_HBV_sites {
	my $seq = shift;
	my @sites_cor_based_on_1;
	my %tags=(
		"rt166"=> "CGTCCTGGGC",
		"rt169" => "CTTTCGCAAA",
		"rt173" => "CCTATGGGA",
		"rt180" => "TCCGTTTCTC",
		"rt181" => "GTTTCTCTTG",
		"rt184" => "CAGTTT",
		"rt191" => "TTGTTCAGTG",
		"rt194" => "CGTAGG",
		"rt200" => "CCACTGTTTG",
		"rt202" => "TTTGGCTTTC",
		"rt204" => "TTTCAGCTAT",
		"rt207" => "GATGAT",
		"rt213" => "TTGGGGGCCAAG",
		"rt214" => "GGCCAAGTCT",
		"rt215" => "CAAGTCTGTA",
		"rt217" => "TGTACAGCAT",
		"rt221" => "GAGTCCCTT",
		"rt229" => "CAATTTTCTT",
		"rt233" => "TCTCTGGGT",
		"rt236" => "TATACATTTA",
		"rt237" => "ACATTTAAAC",
		"rt238" => "TTTAAACCCT",
		"rt245" => "AAGATGGGGT",
		"rt250" => "CTAAACTTC",
		"rt256" => "CATAATTGGA"
	);

	
	my %sites_cor_based_on_1;
	for my $tag_name (keys %tags) {
		my $tag=$tags{$tag_name};
		if($tag_name=~/rt(\d+)/){
			my $tag_num=$1;
			my $tag=$tags{$tag_name};
			if($seq=~/$tag/){
				my $previous=$`;
				my $matched=$&;
				my $post=$';
				my $seq_before_site="$previous$matched";
				my $coor_site=length($seq_before_site)+1;
				$sites_cor_based_on_1{"rt$tag_num"}=$coor_site;
			}else{
				my @tag_bases=split //, $tag;
				my $match_flag=0;
				for my $base_pos (0..$#tag_bases) {
					my @new_tag_bases=@tag_bases;
					$new_tag_bases[$base_pos]="[ATGC]";
					
					my $new_tag=join("", @new_tag_bases);
					if($seq=~/$new_tag/){
						my $previous=$`;
						my $matched=$&;
						my $post=$';
						my $seq_before_site="$previous$matched";
						my $coor_site=length($seq_before_site)+1;
						$sites_cor_based_on_1{"rt$tag_num"}=$coor_site;
					}
				}
			}

		}
	}
	return %sites_cor_based_on_1;
}



sub getAMPpeaks {
	my $abi_file= shift;

	my $abi= ABI->new("$abi_file");
	my @traceA=$abi->get_trace("A");
	my @traceG=$abi->get_trace("G");
	my @traceC=$abi->get_trace("C");
	my @traceT=$abi->get_trace("T");

	#get peaks for each base
	my ($indexes_ref_A, $trace_values_ref_A)=getpeaks(\@traceA);
	my ($indexes_ref_G, $trace_values_ref_G)=getpeaks(\@traceG);
	my ($indexes_ref_C, $trace_values_ref_C)=getpeaks(\@traceC);
	my ($indexes_ref_T, $trace_values_ref_T)=getpeaks(\@traceT);

	#get window around primary basecall peaks
	my @primarypeaks_0=$abi->get_base_calls;
	my @primarypeaks_1= map {$_+1} @primarypeaks_0;
	my @primarypeaks1_tmp=(0,@primarypeaks_1);
	my $primarypeaks_diff_ref=diff(\@primarypeaks1_tmp);
	my @primarypeaks_diff = @{$primarypeaks_diff_ref};

	# make starts;
	my @starts=map{($primarypeaks_1[$_]) -0.5*($primarypeaks_diff[$_])} 0..$#primarypeaks_1;

	# make stops;
	my @primarypeaks_1_tmp2 = @primarypeaks_1;
	pop @primarypeaks_1_tmp2;
	my @primarypeaks_diff_tmp1 = @primarypeaks_diff;
	shift @primarypeaks_diff_tmp1;

	my @part1 = map {$primarypeaks_1_tmp2[$_] + 0.5* $primarypeaks_diff_tmp1[$_]} 0..$#primarypeaks_1_tmp2;
	my @primarypeaks_1_tmp3 = @primarypeaks_1;
	my @primarypeaks_diff_tmp2 = @primarypeaks_diff;
	my $length_primarypeaks_diff = @primarypeaks_diff_tmp2;
	my $last_diff = pop @primarypeaks_diff_tmp2;
	my $primarypeaks_length_diff = $primarypeaks_1_tmp3[$length_primarypeaks_diff-1];
	my $part2 = $primarypeaks_length_diff + 0.5*$last_diff;
	my @stops=(@part1, $part2);
	# hack for last peak. Just uses distance preceding peak as distance after peak.

	# now get max peak value for each channel in each peak window. If no peak return 0.
	my %tempAmpMatrix;
	my %tempPosMatrix;
	my $ratio = 0.33;
	my $primary="";
	my $secondary="";

	my %degenerates=(
		"A" => ["A"],
		"T" => ["T"],
		"C" => ["C"],
		"G" => ["G"],
		"R" => ["A","G"],
		"Y" => ["C","T"],
		"M" => ["A","C"],
		"K" => ["G","T"],
		"S" => ["G","C"],
		"W" => ["A","T"],
		"H" => ["A","T","C"],
		"B" => ["G","T","C"],
		"V" => ["G","A","C"],
		"D" => ["G","A","T"],
		"N" => ["A","T","C","G"]
	);
	my %degenerate_bases;

	for my $degenerate_base (keys %degenerates) {
		my @bases=@{$degenerates{$degenerate_base}};
		my @bases_sort = sort @bases;
		my $bases=join("", @bases_sort);
		$degenerate_bases{$bases}=$degenerate_base;
	}

	for my $start_index (0..$#starts) {
		my $start_i = $starts[$start_index];
		my $stop_i = $stops[$start_index];
		my ($peak_trace_A_value, $peak_trace_A_index) = peakvalues($indexes_ref_A, $trace_values_ref_A, $start_i, $stop_i);
		my ($peak_trace_C_value, $peak_trace_C_index) = peakvalues($indexes_ref_C, $trace_values_ref_C, $start_i, $stop_i);
		my ($peak_trace_G_value, $peak_trace_G_index) = peakvalues($indexes_ref_G, $trace_values_ref_G, $start_i, $stop_i);
		my ($peak_trace_T_value, $peak_trace_T_index) = peakvalues($indexes_ref_T, $trace_values_ref_T, $start_i, $stop_i);
		
		unless($peak_trace_A_index eq "NA" && $peak_trace_C_index eq "NA" && $peak_trace_G_index eq "NA" && $peak_trace_T_index eq "NA"){
			
			my @signals=($peak_trace_A_value,$peak_trace_C_value,$peak_trace_G_value,$peak_trace_T_value);
			my @positions=($peak_trace_A_index,$peak_trace_C_index,$peak_trace_G_index,$peak_trace_T_index);
			$tempAmpMatrix{$start_index+1}=\@signals;
			$tempPosMatrix{$start_index+1}=\@positions;
			
			my $max_signal = max(\@signals);
			my @signalratios=map {$_/$max_signal} @signals;
			my @Bases=qw(A C G T);
			my @filter_indexes=grep {$signalratios[$_] < $ratio } 0..$#signalratios;
			for my $filter_index (@filter_indexes) {
				$Bases[$filter_index]="NA";
			}
			
			my @sorted_indexes = sort {$signals[$b] <=> $signals[$a]} 0..$#signals;
			#print ("@sorted_indexes\n");
			my @Bases_sorted=@Bases[@sorted_indexes];
			my @positions_sorted=@positions[@sorted_indexes];

			my @na_indexes=grep {$_ eq "NA"} @Bases_sorted;
			my $na_count = @na_indexes;
			my $not_na_count = 4-$na_count;
			if($not_na_count ==4 || $not_na_count ==0){
				$primary.="N";
				$secondary.="N";
			}elsif($not_na_count >1){
				$primary.=$Bases_sorted[0];
				my @Bases2 = @Bases_sorted[1..3];
				my @Bases2_not_NA;
				for my $Base2(@Bases2) {
					unless($Base2 eq "NA"){
						push @Bases2_not_NA, $Base2;
					}
				}
				my @Bases2_not_NA_sort = sort @Bases2_not_NA;
				my $Bases2_not_NA_sort = join("", @Bases2_not_NA_sort);
				my $degenerate_base=$degenerate_bases{$Bases2_not_NA_sort};
				$secondary.=$degenerate_base;
			}else{
				$primary.=$Bases_sorted[0];
				$secondary.=$Bases_sorted[0];
			}

		}
	}
	
	my $amppeaks_long="";

	for my $count_one (sort {$a<=>$b} keys %tempAmpMatrix) {
		my @signals_one = @{$tempAmpMatrix{$count_one}};
		my @sort_signals_one = sort {$b<=>$a} @signals_one;
		my $ratio_2 = $sort_signals_one[1]/$sort_signals_one[0];
		my $sig="";
		if($ratio_2>0.2){
			$sig="TRUE";
		}else{
			$sig="FALSE";
		}
		push @signals_one, $ratio_2;
		push @signals_one, $sig;

		my $signals_one=join("\t", @signals_one);
		$amppeaks_long.="$signals_one\n";
	}
	return $amppeaks_long;
}



sub max {
	my $data_ref = shift;
	my @data = @{$data_ref};
	my $max=$data[0];
	for my $data (@data) {
		if($data>$max){
			$max=$data;
		}
	}
	return $max;
}



sub peakvalues {
	my $indexes_ref = shift;
	my $trace_values_ref = shift;
	my $pstart = shift;
	my $pstop = shift;
	my @indexes=@{$indexes_ref};
	my @trace_values = @{$trace_values_ref};

	my @region = grep {$indexes[$_] > $pstart && $indexes[$_]< $pstop} 0..$#indexes;

	if (@region == 0) {
		return (0,"NA");
	}else{
		my @indexes_tmp=@indexes[@region];
		#print ("@indexes_tmp");
		my $max_index = $region[0];
		for my $index_tmp (@region) {
			if($trace_values[$index_tmp] > $trace_values[$max_index]){
				$max_index = $index_tmp;
			}
		}
		return($trace_values[$max_index], $indexes[$max_index]);
	}
}



sub getpeaks { 
	my $trace_ref = shift;
	my ($values_ref, $lengths_ref) = rle($trace_ref);
	my @values = @{$values_ref};
	my @lengths = @{$lengths_ref};
	unshift @values, "-inf";
	push @values, "-inf";
	my $diff1_ref = diff(\@values);
	my @sign_diff1 = sign(@{$diff1_ref});
	my $diff2_ref = diff(\@sign_diff1);
	my @diff3 = map {$_ == -2 ? "TRUE" : "FALSE"}  @{$diff2_ref};
	my @rep=map {($diff3[$_]) x $lengths[$_]} 0..$#diff3;
	my @indexes = map {$_+1} grep { $rep[$_] eq "TRUE" } 0..$#rep ;
	my @indexes1= map {$_-1} @indexes;
	my @trace_values=@{$trace_ref}[@indexes1];
	return (\@indexes, \@trace_values);
}



sub diff {
	my $data_ref = shift;
	my @data=@{$data_ref};
	my @diff;
	for my $index (1..$#data) {
		my $current_value = $data[$index];
		my $previous_value = $data[$index-1];
		my $diff = $current_value-$previous_value;
		push @diff, $diff;
	}
	return \@diff;
}


sub rle {
	my $data_ref = shift;
	my @data = @{$data_ref};
	my @values;
	my @lengths;
	my $first_value = shift @data;
	push @values, $first_value;
	push @lengths, 1;
	for my $value (@data) {
		my $previous_value = $values[$#values];
		my $previous_length = $lengths[$#lengths];
		unless($previous_value eq $value){
			push @values, $value;
			push @lengths, 1;
		}
		else{
			$lengths[$#lengths]++;
		}
	}
	return (\@values, \@lengths);
}



sub sign
{
	return wantarray? map{($_ < 0)? -1: (($_ > 0)? 1: 0)} @_:
		($_[0] < 0)? -1: (($_[0] > 0)? 1: 0);
}


sub plotchroma {

	# === CONFIGURATION ===
	my $abi_file = shift;
	my $prefix;
	if($abi_file=~/ABI\/(.+)\.ab1/){
		$prefix=$1;
	}
	my $output_image="$prefix.chromatogram.png";

	# === LOAD ABI FILE ===
	my $abi = ABI->new($abi_file);
	die "Could not open ABI file '$abi_file'\n" unless $abi;

	# === FETCH TRACE DATA ===
	my @A_trace = $abi->get_trace("A");
	my @C_trace = $abi->get_trace("C");
	my @G_trace = $abi->get_trace("G");
	my @T_trace = $abi->get_trace("T");

	my $trace_length = $abi->get_trace_length();    # X-axis range
	my $sequence     = $abi->get_sequence();        # Full base-called sequence as a string
	my @positions    = $abi->get_base_calls();      # Array of peak positions

	my @bases = split('', $sequence);               # Convert sequence string to array

	# === SETUP IMAGE ===
	my $width  = $trace_length;
	my $height = 450;

	my $img = GD::Image->new($width, $height);
	my $black  = $img->colorAllocate(0, 0, 0);
	my $white  = $img->colorAllocate(255, 255, 255);
	my $green  = $img->colorAllocate(0, 255, 0);
	my $blue   = $img->colorAllocate(0, 0, 255);
	my $red    = $img->colorAllocate(255, 0, 0);
	my $yellow = $img->colorAllocate(255, 255, 0);

	$img->fill(0, 0, $white);  # Black background
	# === NORMALIZE TRACE HEIGHT ===
	my $max_val = 1;
	foreach my $v (@A_trace, @C_trace, @G_trace, @T_trace) {
		$max_val = $v if $v > $max_val;
	}

	# === DRAW TRACE LINES ===
	for my $x (1 .. $trace_length - 1) {
		my $scale = 300;  # Trace vertical scale
		my $y_offset = $height - 60;

		$img->line($x - 1, $y_offset - int($A_trace[$x - 1] / $max_val * $scale),
				   $x,     $y_offset - int($A_trace[$x]     / $max_val * $scale), $green);

		$img->line($x - 1, $y_offset - int($C_trace[$x - 1] / $max_val * $scale),
				   $x,     $y_offset - int($C_trace[$x]     / $max_val * $scale), $blue);

		$img->line($x - 1, $y_offset - int($G_trace[$x - 1] / $max_val * $scale),
				   $x,     $y_offset - int($G_trace[$x]     / $max_val * $scale), $black);

		$img->line($x - 1, $y_offset - int($T_trace[$x - 1] / $max_val * $scale),
				   $x,     $y_offset - int($T_trace[$x]     / $max_val * $scale), $red);
	}

	# === DRAW BASE LETTERS ===
	my $font = gdSmallFont;
	my $text_y = $height - 35;

	for (my $i = 0; $i < @positions && $i < @bases; $i++) {
		my $x = $positions[$i];
		my $base = $bases[$i];

		my $color = $white;
		$color = $green  if $base eq 'A';
		$color = $blue   if $base eq 'C';
		$color = $black if $base eq 'G';
		$color = $red    if $base eq 'T';

		# Center text slightly to left
		$img->string($font, $x - 3, $text_y, $base, $color);
	}

	# === WRITE IMAGE FILE ===
	open my $out, '>', "chromatogram/$output_image" or die "Cannot open $output_image: $!\n";
	binmode $out;
	print $out $img->png;
	close $out;

	print "Chromatogram image for $prefix is saved to $output_image\n";
}
