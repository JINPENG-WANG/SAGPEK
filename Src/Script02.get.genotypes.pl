#!/usr/bin/perl -w
use strict;
use IO::File;

my $type=shift @ARGV;
my $orientation = "f";


my %custom;
if($type=~/custom/i){
	my $tag_file=shift @ARGV;
	if (@ARGV>0) {
		$orientation = shift @ARGV;
		print "   Orientation of Sanger sequencing is Reverse!\n";

	}
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

my @ampfiles=<Amp/*peakAmp.txt>;
if(@ampfiles <1){
	die(" Error: No intermediate peakAmp files found! Please make sure ABI-format files are listed under the ABI directory!\n");
}

my $fh_geno_out;

if($type eq "TEST"){
	$fh_geno_out=IO::File->new(">Genotype/TEST.genotype.txt");
	$fh_geno_out->print("SampleID");
	$fh_geno_out->print("\tFirst_Site_Genotype\tFirst_site_AA\tSecond_Site_Genotype\tSecond_site_AA\tThird_Site_Genotype\tThird_Site_AA\n");
}

if($type eq "cattle_CN"){
	$fh_geno_out= IO::File->new(">Genotype/cattle_CN.genotype.txt");
	$fh_geno_out->print("SampleID\tCN_genotype\tAmino Acids\n");
}

if($type eq "cattle_DUMPS"){
	$fh_geno_out=IO::File->new(">Genotype/cattle_DUMPS.genotype.txt");
	$fh_geno_out->print("SampleID\tDUMPS_genotype\tAmino Acids\n");
}

if($type eq "HBV"){
	$fh_geno_out=IO::File->new(">Genotype/HBV.mutants.txt");
	$fh_geno_out->print("SampleID");
	my @rt_nums=(166, 169, 173, 180, 181, 184, 191, 194, 200, 202, 204, 207, 213, 214, 215, 217, 221, 229, 233, 236, 237, 238, 245, 250, 256);
	for my $rt_num(@rt_nums) {
		$fh_geno_out->print("\trt$rt_num");
	}
	$fh_geno_out->print("\n");
}

if($type eq "PAH"){
	$fh_geno_out = IO::File->new(">Genotype/Newborn_PAH.mutants.txt");
	$fh_geno_out->print("SampleID");
	my @AA_nums = (53, 107, 111, 171, 223, 225, 243, 304, 356, 388, 392, 399, 408, 413);
	for my $AA_num (@AA_nums) {
		$fh_geno_out->print("\tAA_${AA_num}_condons\tAA_${AA_num}_Amino_Acids");
	}
	$fh_geno_out->print("\n");
}

if($type=~/custom/i){
	$fh_geno_out = IO::File->new(">Genotype/Custom.mutants.txt");
	$fh_geno_out ->print("SampleID");
	for my $custom_tag_count (sort {$a<=>$b} keys %custom) {
		my $custom_tag_name = $custom{$custom_tag_count}{tag_name};
		$fh_geno_out->print("\t$custom_tag_name-condons\t$custom_tag_name-AA");
	}
	$fh_geno_out->print("\n");
}
print " #SAGPEK: Step 2: Performing genotyping.\n"; 
for my $ampfile (@ampfiles) {
	my $fh_in = IO::File->new("$ampfile",'r');
	my $prefix;
	if($ampfile=~/^Amp\/(.+)\.ab1\.peakAmp\.txt/){
		$prefix=$1;
		print "   Processing $prefix......";
	}
	$fh_geno_out->print("$prefix");
	my $line_count=0;
	my $primary_seq;
	my $secondary_seq;
	my @sigs;
	my %complementary = (
		"A" => "T",
		"T" => "A",
		"C" => "G",
		"G" => "C"
	);

	while(<$fh_in>){
		chomp;
		my $line = $_;
		$line_count++;
		my %data;
		if($line_count>1){
			my ($A,$C,$G,$T,$ratio,$sig)=split /\s+/, $line;
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
	}

	### TMP -1 
	#my $fh_out111=IO::File->new(">Genotype/$prefix.seqs.txt");
	#$fh_out111->print(">primary_seq\n$primary_seq\n");
	#$fh_out111->print(">secondary_seq\n$secondary_seq\n");

	if($orientation eq "r"){
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
	
	#####  TMP-2
	#$fh_out111->print(">r_primary\n$primary_seq\n");
	#$fh_out111->print(">r_secondary\n$secondary_seq\n");
	######
	
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
	if($type eq "TEST"){	
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
			$fh_geno_out->print("\t$condon_all_possible_string\t$AA_all_possible_string");
		}
		$fh_geno_out->print("\n");
	}


	#####...#####
	# For HBV.
	if($type eq "HBV"){
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
	if($type =~/cattle/){

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
	if($type =~/PAH/){
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
	if($type =~/custom/i){
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
				}
			}
			$fh_geno_out->print("\t$condon_all_possible_string\t$AA_all_possible_string");
		}
		$fh_geno_out->print("\n");
	}

	print "Done!\n";
}

print " #SAGPEK: Whole program completed!\n";


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
