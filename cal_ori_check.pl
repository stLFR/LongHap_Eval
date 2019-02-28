#! /usr/bin/perl
use strict;

my $var = $ARGV[0];
my $num_ori;
my $num_ori_snp;
my $num_ori_indel;
my %hash;
open VAR, $var or die $!;
while(<VAR>){
	chomp;
	next if /^#/;
	my ($ref, $pos, $ref_base, $alt_base) = (split /\s+/, $_)[0,1,3,4];
	next if /LowQual/;
	next unless ($ref =~ /chr\d+/ or $ref eq "chrX");
	my $info = (split /\s+/, $_)[-1];
	my $depth = (split /:/, $info)[2];
	if ($info=~/1\/2/ or $info=~ /2\/1/ or $info=~/1\|2/ or $info=~/2\|1/){ 
				$num_ori++;
				$hash{$ref}++;
				$alt_base =~ /(\w+),(\w+)/;
				my ($alt_1, $alt_2) = ($1, $2);
				my $alt_1_len = length $alt_1;
				my $alt_2_len = length $alt_2;
				my $ref_len = length $ref_base;
				unless ($alt_1_len eq $ref_len and $alt_2_len eq $ref_len){
					$num_ori_indel++;
				}else{
					$num_ori_snp++;
				}
	}elsif ($info=~/0\/1/ or $info=~/0\|1/ or $info=~/1\|0/){
				$num_ori++;
				$hash{$ref}++;
				my $alt_len = length $alt_base;
				my $ref_len = length $ref_base;
				unless ($alt_len eq $ref_len){
					$num_ori_indel++;
				}else{
					$num_ori_snp++;
				}
	}
}
close VAR;
print $num_ori, "\n", $num_ori_snp, "\n", $num_ori_indel, "\n";

my @files = `find $ARGV[1] -name '*.log'`;
chomp @files;

my $total;
my $total_snp;
my $total_indel;
my %hash2;
my %snp;
my %indel;
foreach my $file (@files){
	open IN, $file or die $!;
	$/ = "break";
	<IN>;
	while(<IN>){
		my ($line1, $line2) = (split /\n/, $_)[1,2];
		my @tmp1 = (split /;/, $line1);
		my @tmp2 = (split /;/, $line2);
		my $num = @tmp1;
		next unless $num >1;
		my %check;
		for (my $i=0; $i<@tmp1; $i++){
			$tmp1[$i]=~/(.+)_(\d+)_(\w+)/;
			my ($chr, $pos, $alle1) = ($1, $2, $3);
			$tmp2[$i]=~/(.+)_(\d+)_(\w+)/;
			my $alle2 = $3;
			#next unless $chr eq "chr1";
			if ($alle1 =~ /[I|D]/ or $alle2 =~ /[I|D]/){
				$indel{$chr}{$pos} =1;
				$total++;
				$total_indel++;
				$hash2{$chr}++;
			}elsif($alle1 =~ /[A|T|G|C]/ and $alle2 =~ /[A|T|G|C]/){
				$snp{$chr}{$pos} =1;
				$total++;
				$total_snp++;
				$hash2{$chr}++;
			}else{
				die "wrong format in $_\n";
			}
			my @a = split /_/, $tmp2[$i];
                        my $b = @a;
                        die "wrong format in $_\n" if $b ne 3;
		}
	}
	close IN;
}
my $phased_snp;
foreach my $chr (keys %snp){
	foreach my $pos (keys %{$snp{$chr}}){
		$phased_snp++;	
	}
}
my $phased_indel;
foreach my $chr (keys %indel){
        foreach my $pos (keys %{$indel{$chr}}){
                $phased_indel++;
        }
}
print "phased_snp:$phased_snp\tphased_indel:$phased_indel\n";

print $total, "\n", $total_snp, "\n", $total_indel, "\n";
print $total/$num_ori, "\n";
print $total_snp/$num_ori_snp, "\n";
print $total_indel/$num_ori_indel, "\n";

foreach my $chr (keys %hash2){
	print $chr, "\t", $hash2{$chr}, "\t", $hash{$chr}, "\t", $hash2{$chr}/$hash{$chr}, "\n";
}
