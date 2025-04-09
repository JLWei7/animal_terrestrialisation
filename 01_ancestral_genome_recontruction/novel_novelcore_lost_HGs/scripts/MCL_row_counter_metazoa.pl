#!/usr/bin/perl -w
#Jordi Paps, Oxford 2015
#Beginning Perl for Bioinformatics
use strict;
use Term::ANSIColor;

#Script to parse the output of MCL, checking for presence/absence of each species in each group of homology
#Each line of the MCL output is an homology group, each line will be checked for the presence of species labels/tags

#Program USAGE
my $USAGE = "Script to parse the output of MCL, checking for presence/absence of each species in each group of homology
USAGE: $0 MCL_output\n";
unless (@ARGV) {
	print color ("red"), $USAGE, color ("reset") and exit;
}

#Check if the file exists, open it
my $filename = $ARGV[0];
unless (check_file ($filename)) {
	open(FILE, $filename) or die "Could not open file $filename\n";
}

#Open output file
my $outfile = $filename;
$outfile =~ s/...(.*)\....$/$1/g;
$outfile = "02_".$outfile."_gene_numbers_parsed.out";
print "Outfile name: $outfile\n";
unless ( open(OUTPUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}

#Count each taxa ocurrence in each line,
my @taxa = array_taxa ();
print "Taxa: @taxa\n";
my $line_number = 0;
my @output = ();

while (<FILE>) {
	$line_number++;
	@output = ();
	push (@output, $line_number);
	foreach my $taxon (@taxa){
		my $count = 0;
		while ($_ =~ /$taxon/g) {
			$count++;
		}
		push (@output, $count);
	}
	print OUTPUT "@output\n"; #Print results to output file, one line at a time.
}
print "\nResults printed in file $outfile\n";

#Close file handles
close(OUTPUT);
close FILE;

#End of program
print color ("red"), "\nend\n\n", color ("reset");
exit;


################################################################################
#                                Subroutines                                   #
################################################################################

#CHECK_FILE: checks file properties, checks if file exists, is flat, is empty, or can be opened
#-e exists, -f flat file, -s empty
sub check_file {
    my ($filename) = @_;
    unless (-e $filename || (undef $filename)) {print "File does NOT exists\n" and exit;}
    unless (-f $filename) {print "File is NOT a flat file\n" and exit;}
    unless (-s $filename) {print "File is empty\n" and exit;}
    unless (open (FH, $filename)) {print "Can not be opened\n" and exit;}
    close FH;
    return;
}

#ARRAY_TAXA: just an array definition to not clutter the main program.
sub array_taxa {
	my (@taxa) = qw/Cfra Cowc Saro Aque Emue Mlei Tadh Aequ Nvec Hvul Hmia Bcal Aric Aste Dcar Rsoc Rsor Rota Mlig Pcro Smed Lana Paus Bner Pesc Dgyr Lluy Pech Ctel Eand Mvul Hrob Hman Hmed Llon Hhan Obim Spha Dpol Mmer Mcal Myes Cgig Lgig Ppel Hrub Gmag Gaeg Mcor Acal Echl Bgla Aful Cuni Oida Pcau Tuco Psam Hali Ppac Dcor Cele Rhab Anan Pred Pdav Hexe Rvar Ebro Nstr Lpol Popi Dpte Dsil Cscu Ptep Gmae Hhol Smar Rimm Dste Eaff Tcal Ppol Pchi Hame Pcla Cqua Anas Hazt Dmag Dgal Dpul Dapu Tqin Caug Iele Znev Tcas Dmel Skow Anja Ajap Spur Arub Blan Cint Ebur Cmil Leri Cpla Drer Sorb Bpec Pmag Pmod Onil Olat Trub Lcha Pann Rbiv Gser Muni Pwal Bbom Xlae Xtro Sbom Bbuf Npar Rtem Emac Acar Ctig Psin Deco Care Cmyd Cser Pmeg Cpic Tcar Amis Asin Cpor Ggan Ggal Tgut Oana Mdom Hsap Clup Ttru/;
	return @taxa;
}
