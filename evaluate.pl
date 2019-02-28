#! /usr/bin/perl
use strict;

my @vcfs = `find $ARGV[0] -name "*.vcf"`;
chomp @vcfs;

my %hash;
my (%hap0, %hap1, %check);
foreach my $vcf (@vcfs){
	open IN, $vcf or die $!;
	while(<IN>){
		chomp;
		next if /^#/;	
		my ($chr, $pos, $ref, $alt) = (split /\s+/, $_)[0,1,3,4];
		my ($FORMAT_DES, $FORMAT) = (split /\s+/, $_)[-2,-1]; #NOTE: the last column could be 1|1:614:0,264:99,361:400:1/1:.:PATMAT or 0|1:76485_G_T:631:125,112:169,167:1007
		my @DES = split /:/, $FORMAT_DES;
		my ($GT, $IGT, $set);
		for (my $i=0; $i<@DES; $i++){
			my $tmp = $DES[$i];
			my $tmp_value = (split /:/, $FORMAT)[$i];
			if ($tmp eq "GT"){
				$GT = $tmp_value;
			}elsif ($tmp eq "IGT"){
				$IGT = $tmp_value;
			}elsif ($tmp eq "PS"){
				$set = $tmp_value;
			}
		}
		$hash{$set}=1;
		warn "1/0 in $_!\n" if $IGT eq "1/0";
		$pos = "$chr:$pos";
#		if ($IGT eq "0/1" or ($GT eq $IGT and $IGT eq "0|1") or ($GT eq $IGT and $IGT eq "1|0")){
		if ($GT eq "0/1" or $GT eq "0|1" or $GT eq "1|0"){
			next if ($IGT eq "0|1" and $IGT ne $GT) or ($IGT eq "1|0" and $IGT ne $GT);
			my $nt;
			($nt, $alt) = alle($alt, $ref);
			$ref=~ /^([A|T|G|C])/;
                	$ref= $1;
			if ($GT eq "0|1"){
				$hap0{$set}{$pos} = $ref;
				$hap1{$set}{$pos} = $nt;
				$check{$pos} =$set;
			}elsif ($GT eq "1|0"){
				$hap0{$set}{$pos} = $nt;
				$hap1{$set}{$pos} = $ref;
				$check{$pos} =$set;
			}else{
				next;				
			}
#		}elsif ($IGT eq "1/2" or $IGT eq "2/1" or ($GT eq $IGT and $IGT eq "1|2") or ($GT eq $IGT and $IGT eq "2|1")){
		}elsif ($GT eq "1/2" or $GT eq "1|2" or $GT eq "2|1"){
			next if ($IGT eq "1|2" and $IGT ne $GT) or ($IGT eq "2|1" and $IGT ne $GT);
			$alt =~ /(\w+),(\w+)/;
                	my ($alle_1, $alle_2) = ($1, $2);
			my ($nt_1, $nt_2);
                	($nt_1, $alle_1) = alle ($alle_1, $ref);
			($nt_2, $alle_2) = alle ($alle_2, $ref);
#			next unless ((length $alt) ==1 and (length $alt) ==2);
			if ($GT eq "1|2"){
				$hap0{$set}{$pos} = $nt_1;
                                $hap1{$set}{$pos} = $nt_2;
				$check{$pos} =$set;
			}elsif ($GT eq "2|1"){
				$hap0{$set}{$pos} = $nt_2;
                                $hap1{$set}{$pos} = $nt_1;
				$check{$pos} =$set;
			}
		}
	}
	close IN;
}

sub alle {
        my ($alt, $ref) = @_;
        my $nt_of_alt;
        my $alt_len = length $alt;
        my $ref_len = length $ref;
        if ($alt_len < $ref_len){     #del
                my $indel_len = $ref_len - $alt_len;
                $nt_of_alt = substr ($ref, $alt_len, $indel_len);
                $nt_of_alt = "D$nt_of_alt";
                $alt = "D$indel_len";
        }elsif ($alt_len > $ref_len){ #ins
                my $indel_len = $alt_len - $ref_len;
                $nt_of_alt = substr ($alt, $ref_len, $indel_len);
                $nt_of_alt = "I$nt_of_alt";
                $alt = "I$indel_len";
        }else{                        #snp
                $alt =~ /^([A|T|G|C])/; #GA G,AA
#                $nt_of_alt = $alt;
                $alt = $1;
		$nt_of_alt = $alt;
        }
        return ($nt_of_alt, $alt);
}


my @files = `find $ARGV[1] -name "*.log"`;
chomp @files;

$/ = "break";
my ($total, $short_snp, $short_indel, $short_alle, $long_comb, $long_snp, $long_indel,$fp_snp, $fp_indel) = (0,0,0,0,0,0,0,0,0);
foreach my $file (@files){
	open IN, $file or die $!;
	<IN>;
	while(<IN>){
		chomp;
		my ($out0, $out1);
		my ($num_of_0, $num_of_1);
		my ($line0, $line1) = (split /\n/,$_)[1,2];
		my @block0 = split /;/, $line0;
		my @block1 = split /;/, $line1;
		my %hash0;
		my %hash1;
		for (my $i=0; $i<@block0; $i++){
			$block0[$i] =~ /(\w+)_(\w+)_(\w+)/;
			my ($chr, $pos, $base0) = ($1, $2, $3);
			next unless ($check{"$chr:$pos"});
			$total++;
			$block1[$i] =~ /(\w+)_(\w+)_(\w+)/;
			my $base1 = $3;
			$hash0{$chr}{$pos}= $base0;
			$hash1{$chr}{$pos}= $base1;
		}
		my (%out0, %out1, %num0, %num1);
		foreach my $chr (keys %hash0){
			foreach my $pos (sort {$a<=>$b} keys %{$hash0{$chr}}){		
				my $base0 = $hash0{$chr}{$pos};
				my $base1 = $hash1{$chr}{$pos};
				$pos = "$chr:$pos";
				my $set = $check{$pos};
				if ($base0 eq $hap0{$set}{$pos}){
					if (length $base0 ==1){
						$out0{$set} .= "0"; 
					}else{
						$out0{$set} .= "2";
					}
					$num0{$set} ++;				
				}elsif ($base0 eq $hap1{$set}{$pos}){
					if (length $base0 ==1){
						$out0{$set} .= "1";
					}else{
						$out0{$set} .= "3";
					}
					$num1{$set} ++;
				}else{
					$out0{$set} .= "4";
					if ((length $base0) ==1 and (length $base1)==1){
						$fp_snp++;
					}else{
						$fp_indel++;
					#	warn $chr, "\t", $pos, "\t", $base0, "\t", $base1, "\n";
					}
				}
				if ($base1 eq $hap1{$set}{$pos}){
					if (length $base1 ==1){
						$out1{$set} .= "1";
					}else{
						$out1{$set} .= "3";
					}
				}elsif ($base1 eq $hap0{$set}{$pos}){
					if (length $base1 ==1){
						$out1{$set} .= "0";
					}else{
						$out1{$set} .= "2";
					}
				}else{
					$out1{$set} .= "4";
					if ((length $base0) ==1 and (length $base1)==1){
                                	        $fp_snp++;
	                                }else{
        	                                $fp_indel++;
	#					warn $chr, "\t", $pos, "\t", $base0, "\t", $base1, "\n";
        	                        }
				}
			}
		}
		foreach my $set (keys %out1){
			next unless $out1{$set} =~ /\w+/;
			my $out0 = $out0{$set};
			my $out1 = $out1{$set};
			my @info; my @info_another;
			if ($num0{$set} >= $num1{$set}){
				$out0 =~ s/2/0/g;
				$out1 =~ s/3/1/g;
				@info = split /0+/, $out0;
				@info_another = split /1+/, $out1;
			}else{
				$out0 =~ s/3/1/g;
                	        $out1 =~ s/2/0/g;
				@info = split /0+/, $out1;
				@info_another = split /1+/, $out0;
			}
			for (my $i=0; $i<@info; $i++){
				if ($info[$i] eq "1"){
					$short_snp++;
				}elsif ($info[$i] eq "3"){
					$short_indel++;
				}elsif($info[$i] eq "4"){
					$short_alle++;	
				}elsif ($info[$i] =~ /\d+/){
					if ($info[$i] =~ /3/){
						if ($info[$i] =~ /1/){
							$long_comb++;
						}else{
							$long_indel++;
						}
					}else{
						$long_snp++;
					}
				}
			}
			for (my $i=0; $i<@info_another; $i++){
				if ($info_another[$i] eq "4"){
					$short_alle++;
				}
			}	
		}
#		print $short, "\t", $long, "\n";
	}
	close IN;
}
$/ = "\n";
print $total, "\t", $short_snp, "\t", $short_indel, "\t",$short_alle, "\t", $long_comb, "\t", $long_snp, "\t", $long_indel, "\t",$fp_snp, "\t", $fp_indel, "\n";
print "short switch error:\t", ($short_snp+$short_indel+$short_alle)/$total, "\n"; 
print "short switch error by snp:\t", $short_snp/$total, "\n";
print "short switch error by indel:\t", $short_indel/$total, "\n";
print "short switch error by incorrect allele:\t",$short_alle/$total, "\n";
print "long switch error:\t",($long_comb+$long_snp+$long_indel)/$total, "\n";
print "long switch error by snp:\t",$long_snp/$total, "\n";
print "long switch error by indel:\t",$long_indel/$total, "\n";
print "long switch error mixing snp and indel:\t",$long_comb/$total, "\n";

foreach my $set(keys %hash){
#	warn $set, "\n";
}
