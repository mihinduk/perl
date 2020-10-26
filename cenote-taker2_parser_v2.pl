#!/usr/bin/env perl
# cenote-taker2_parser_v2.pl
# 2020_07_01
use warnings;
use strict;

# This script takes the output from Cenote-taker2 and parses it to generate an outfile with taxonomy for R

# $ARGV[0] = output from Cenote-taker2
# $ARGV[1] = Cenote-taker2 mode: DNA (default) or RNA

# Make $ARGV[1] lc
my $na = lc($ARGV[1]); #print "\$na = $na\n";

my $outfile = $ARGV[0];
$outfile =~ s/\.tsv/_clean_tax.txt/; #print "\$outfile = $outfile\n";

my %non_euk_family=();
my %non_euk = ();
my ($unneeded, $to_be_split, $GI, $ref, $to_be_split_again, $Species_to_find, $ne_key, $ne_val, $protein);
my $x = 0;

# Read the taxonomy files into a hash
# non_euk
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_non_euk_fin_fib.txt") or die "Couldn't find the non_euk file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $key = "$tax_id;$Kingdom;$Phylum;$Class;$Order;$Family;$Genus;$Species";
	$non_euk {$key} = $tax_name; #print "\$tax_name $tax_name linked to $key\n";
}
close(IN);

# Viruses
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_viruses_fin_fib.txt") or die "Couldn't find the viral file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $v_upstream = "$Kingdom".";$Phylum".";$Class".";$Order";
	$non_euk_family{$Family} = $v_upstream;
}
close(IN);

# Archaea
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_archaea_fin_fib.txt") or die "Couldn't find the viral file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $a_upstream = "$Kingdom".";$Phylum".";$Class".";$Order";
	$non_euk_family{$Family} = $a_upstream; #print "$Family = $a_upstream\n";
}
close(IN);

# Bacteria
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_bacteria_fin_fib.txt") or die "Couldn't find the viral file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $b_upstream = "$Kingdom".";$Phylum".";$Class".";$Order";
	$non_euk_family{$Family} = $b_upstream;
}
close(IN);

# Parse the DNA file
if ($na eq "dna") {
	open(IN, $ARGV[0]) or die "Couldn't find the cenote-taker2 file $!";
	open (OUT, ">$outfile") or die "Couldn't create a file for your clean taxonomy $!";
	print OUT "Isolation source\tCompleteness\tCenote-taker contig name\toriginal contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\n";
	while (my $line = <IN>) {
		next if ($line =~/^Isolation source/);
		next if ($line =~/^\t/);
	  chomp ($line);
		my ($IS, $Completeness, $CTcn, $ocn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTP, $BLASTN) = split(/\t/, $line, 11);
		$IS=~ s/^\s+|\s+$//g; $Completeness=~ s/^\s+|\s+$//g; $CTcn=~ s/^\s+|\s+$//g; $ocn=~ s/^\s+|\s+$//g; $Length=~ s/^\s+|\s+$//g;
		$Element=~ s/^\s+|\s+$//g; $Topology=~ s/^\s+|\s+$//g; $CVD=~ s/^\s+|\s+$//g; $ORFc=~ s/^\s+|\s+$//g; $BLASTP=~ s/^\s+|\s+$//g;
		$BLASTN=~ s/^\s+|\s+$//g;
		print OUT "$IS\t$Completeness\t$CTcn\t$ocn\t$Length\t$Element\t$Topology\t$CVD\t$ORFc\t$BLASTN\t";
		if ($BLASTP !~ /;/) {
			#print "Here's something odd\n"; print "\$BLASTP = $BLASTP\n";
			my ($possible_family, $extra) = split (/ /, $Element, 2); #print "\$possible_family = $possible_family\n";
			if (exists($non_euk_family{$possible_family})) {
				my $family_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$possible_family";
				my $genus_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$possible_family"."_Genus_undefined";
				print OUT "GI_undefined\tref_undefined\ttaxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\t$family_out\t$genus_out\t$Element\tno protein\n";
			} else {
			print OUT "GI_undefined\tref_undefined\ttaxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_undefined\t$Element\tno protein\n";
			}
		} else {
			($unneeded, $to_be_split) = split(/gi\|/, $BLASTP, 2); #print "\$to_be_split = $to_be_split\n";
			($GI, $unneeded, $ref, $to_be_split_again) = split(/\|/, $to_be_split, 4); #print "\$to_be_split_again = $to_be_split_again\t\$GI = $GI\t\$ref = $ref\n";
			print OUT "$GI\t$ref\t";
			($protein, $to_be_split) = split(/\[/, $to_be_split_again, 2);  #print "\$to_be_split = $to_be_split\n";
			$protein=~ s/^\s+|\s+$//g;
			($Species_to_find, $unneeded) = split(/\]/, $to_be_split, 2);  #print "\$Species_to_find = $Species_to_find\n";
			my @non_euk_matches = grep /.*$Species_to_find*/, values %non_euk; #print @matches; print "\n";
			if (scalar(@non_euk_matches) >= 1) {
				my $non_euk_length = scalar(@non_euk_matches); #print "\$non_euk_length = $non_euk_length\n";
				if ($non_euk_length == 1) {
					$ne_val = shift(@non_euk_matches); #print "\$ne_val = $ne_val\n";
					($ne_key) = grep { $non_euk{$_} eq $ne_val } keys %non_euk; #print "\$ne_key = $ne_key\n";
					my ($tax_id, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $ne_key, 8); #print "\$tax_id = $tax_id\n";
					print OUT "$tax_id\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t$protein\n";
				} else {
					#print "Need to figure these out contigs\n"; print @non_euk_matches; print "\n";
					$ne_val = shift(@non_euk_matches); #print "\$ne_val = $ne_val\n";
					($ne_key) = grep { $non_euk{$_} eq $ne_val } keys %non_euk; #print "\$ne_key = $ne_key\n";
					my ($tax_id, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $ne_key, 8); #print "\$tax_id = $tax_id\n";
					my $Species_fin = "$Species" . "_multi"; #print "\$Species_fin = $Species_fin\n";
					print OUT "$tax_id\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species_fin\t$protein\n";
				}
			} elsif (scalar(@non_euk_matches) == 0) {
					#print "Nothing matches $Species_to_find\n";
					print OUT "taxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_undefined\t$Species_to_find\t$protein\n";
				}
			}
		}
	} elsif ($na eq "rna") {
		open(IN, $ARGV[0]) or die "Couldn't find the cenote-taker2 file $!";
		open (OUT, ">$outfile") or die "Couldn't create a file for your clean taxonomy $!";
		print OUT "Isolation source\tCompleteness\tCenote-taker contig name\toriginal contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\n";
		while (my $line = <IN>) {
			next if ($line =~/^Isolation source/);
			next if ($line =~/^\t/);
			chomp ($line);
			my ($IS, $Completeness, $CTcn, $ocn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTP, $BLASTN) = split(/\t/, $line, 11);
			$IS=~ s/^\s+|\s+$//g; $Completeness=~ s/^\s+|\s+$//g; $CTcn=~ s/^\s+|\s+$//g; $ocn=~ s/^\s+|\s+$//g; $Length=~ s/^\s+|\s+$//g;
			$Element=~ s/^\s+|\s+$//g; $Topology=~ s/^\s+|\s+$//g; $CVD=~ s/^\s+|\s+$//g; $ORFc=~ s/^\s+|\s+$//g; $BLASTP=~ s/^\s+|\s+$//g;
			$BLASTN=~ s/^\s+|\s+$//g;
			my ($not_wanted, $true_contig) = split(/\@/, $ocn, 2); #print "\$true_contig = $true_contig\n";
			print OUT "$IS\t$Completeness\t$CTcn\t$true_contig\t$Length\t$Element\t$Topology\t$CVD\t$ORFc\t$BLASTN\t";
			if ($BLASTP !~ /;/) {
				#print "Here's something odd\n"; print "\$BLASTP = $BLASTP\n";
				my ($possible_family, $extra) = split (/ /, $Element, 2); #print "\$possible_family = $possible_family\n";
				if (exists($non_euk_family{$possible_family})) {
					my $family_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$possible_family";
					my $genus_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$possible_family"."_Genus_undefined";
					print OUT "GI_undefined\tref_undefined\ttaxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\t$family_out\t$genus_out\t$Element\tno protein\n";
				} else {
				print OUT "GI_undefined\tref_undefined\ttaxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_undefined\t$Element\tno protein\n";
				}
			} else {
				($unneeded, $to_be_split) = split(/gi\|/, $BLASTP, 2); #print "\$to_be_split = $to_be_split\n";
				($GI, $unneeded, $ref, $to_be_split_again) = split(/\|/, $to_be_split, 4); #print "\$to_be_split_again = $to_be_split_again\t\$GI = $GI\t\$ref = $ref\n";
				print OUT "$GI\t$ref\t";
				($protein, $to_be_split) = split(/\[/, $to_be_split_again, 2); #print "\$protein = $protein\t\$to_be_split = $to_be_split\n";
				$protein=~ s/^\s+|\s+$//g;
				($Species_to_find, $unneeded) = split(/\]/, $to_be_split, 2);  #print "\$Species_to_find = $Species_to_find\n";
				my @non_euk_matches = grep /.*$Species_to_find*/, values %non_euk; #print @matches; print "\n";
				if (scalar(@non_euk_matches) >= 1) {
					my $non_euk_length = scalar(@non_euk_matches); #print "\$non_euk_length = $non_euk_length\n";
					if ($non_euk_length == 1) {
						$ne_val = shift(@non_euk_matches); #print "\$ne_val = $ne_val\n";
						($ne_key) = grep { $non_euk{$_} eq $ne_val } keys %non_euk; #print "\$ne_key = $ne_key\n";
						my ($tax_id, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $ne_key, 8); #print "\$tax_id = $tax_id\n";
						print OUT "$tax_id\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t$protein\n";
					} else {
						#print "Need to figure these out RNA contigs\n"; print @non_euk_matches; print "\n";
						$ne_val = shift(@non_euk_matches); #print "\$ne_val = $ne_val\n";
						($ne_key) = grep { $non_euk{$_} eq $ne_val } keys %non_euk; #print "\$ne_key = $ne_key\n";
						my ($tax_id, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $ne_key, 8); #print "\$tax_id = $tax_id\n";
						my $Species_fin = "$Species" . "_multi";
						print OUT "$tax_id\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species_fin\t$protein\n";
					}
				} elsif (scalar(@non_euk_matches) == 0) {
						#print "Nothing matches $Species_to_find\n";
						print OUT "taxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_undefined\t$Species_to_find\t$protein\n";
					}
				}
			}
		}


close(IN);
close(OUT);

print "Ingenio maximus, arte rudis. Maximum ingenuity, raw technique. - Ovid\nYour outfile is $outfile\n";
