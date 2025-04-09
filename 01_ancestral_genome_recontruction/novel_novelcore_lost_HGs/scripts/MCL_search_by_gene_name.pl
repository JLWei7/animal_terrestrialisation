#!/usr/bin/perl -w
#Jordi Paps, Oxford 2015
#Beginning Perl for Bioinformatics
use strict;
use AnyDBM_File;
use Fcntl;
use Term::ANSIColor;

#Script to extract groups of homology from a parsed MCL output, based on the taxonomic distribution

#Introduce file names
my $MCL_out = "Input/Orthogroups.txt"; 
my $MCL_columns_parsed = "Input/Orthogroups_gene_numbers_parsed.txt";

#Check if files exists, open them
check_file ($MCL_out);
check_file ($MCL_columns_parsed);

#Create DBM files from the MCL output file for fast lookups.
my %MCL_out;
tie (%MCL_out,"AnyDBM_File", "$MCL_out.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
#%MCL_out = MCL_output_to_hash ($MCL_out);
print "\nMCL out (showing some random lines out of ", scalar (keys %MCL_out), "):\n";
my @MCL_out_keys = keys %MCL_out;
foreach (@MCL_out_keys[0..4]){
	my $first_elements = substr ($MCL_out{$_}, 0, 80);
	print "$_ => $first_elements...\n";
} print "...\n";

#Create DBM files from the MCL parsed file for fast lookups.
my %MCL_columns_parsed;
tie (%MCL_columns_parsed,"AnyDBM_File", "$MCL_columns_parsed.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
#%MCL_columns_parsed = MCL_columns_parsed_to_hash ($MCL_columns_parsed);
print "\nMCL parsed (each column is one species, showing some random lines out of ", scalar (keys %MCL_columns_parsed), "):\n\n";
print "Species: Crei Ppat Smoe Atri Atha Ehux Bnat Rfil Ttra Falb Spun Amac Scer Sarc Cfra Cow_ Mbre Sros Aque Ocar Mley Pbac Tadh Nvec Adig Hmag Gsal Sjap Sman Egra Emul Hmic Avag Cgig Pfuc Lgig Ctel Hrob Tspi Rcul Cele Bmal Smar Isca Smim Mmar Dpul Znev Tcas Dmel Skow Spur Bflo Cint Csav Bsch Odio Drer Xtro Ggal Acar Hsap\n\n";
my @MCL_columns_parsed_keys = keys %MCL_columns_parsed;
foreach (@MCL_columns_parsed_keys[0...4]){
	print "$_ => $MCL_columns_parsed{$_}\n";
} print "...\n";

#Ask user for clade/spp to check genes taxonomic distribution, perform search
my $user_input = '';	#To store the searching criteria fom user
my @arguments = ();	#To store the split of the searching criteria
my %final_searches = ();#To store ALL searches columns and options (present, absent, minus, atleast,...)
my @good_homology_groups_spp_names = ();	#To store the groups of homology fullfilling the search criteria
my @good_homology_groups_spp_names_copy = ();	#Backup
my @good_homology_groups_columns_parsed = ();	#To store the columns of the groups of homology fullfilling the search criteria

OUTER: do {
	print "\nPlease, enter name of the gene to search (\"example\" for some samples, Enter to exit):\n";
	$user_input = <STDIN>;
	chomp $user_input;
	#print "Query: $user_input\n";
	unless ($user_input =~ /^\s*$|^exit|^quit/i) {			#Unless user wants to exit...
		if ($user_input =~ /example|help/i) {
			print_examples();			#Print examples of commands if user requests it
		} else {
			#Here the real search starts, emptying variables for next loop and parsing the user input
			%final_searches = ();			#Empty the hash containing all the search conditions
			@good_homology_groups_spp_names = ();	#Empty the hash containing the results from previous search
			@good_homology_groups_columns_parsed = (); #Empty the hash containing the columns of the groups of homology fullfilling the search criteria
			@arguments = split (" ", $user_input);	#Decompose the user input in different arguments entered by user, each taxa query delimited by a space (" ")
			#print "\@arguments: @arguments\n";
			my @MCL_out_keys = keys %MCL_out;
			my $HG_counter = 0;
			foreach my $homology_group (@MCL_out_keys) {
				my @flags = ();
				#Read the hash containing the MCL output, check the presence of the search items
				foreach my $query (@arguments){ 			#For each line of the MCL output file...
					#print "HG_$homology_group ", substr ($MCL_out{$homology_group},0, 80),"\n";
					#print "Query: $query\n";
					if ($MCL_out{$homology_group} =~ /$query/gi){						
						print "Match in HG $homology_group\n";
						push (@flags, "true");
					}
					else {
						#print "No match in HG $homology_group\n";
						push (@flags, "false");
					}
				}
				unless (grep {$_ eq "false"} @flags) {					
					push (@good_homology_groups_spp_names, "$homology_group\t$MCL_out{$homology_group}\n");
					push (@good_homology_groups_columns_parsed, "$homology_group\t$MCL_columns_parsed{$homology_group}\n");
					$HG_counter++;
				}
			}
			#Report results, ask to save files
			@good_homology_groups_spp_names_copy = @good_homology_groups_spp_names;
			print "\nNumber of groups of homology found: $HG_counter \n";
			if ($HG_counter == 0) { goto OUTER; };
			
			#Save 4 different output files
			print "\nDo you want to see results (gene names, MCL groups, etc) and save in files? Yes/No\n";
			my $save_files = <STDIN>;
			chomp $save_files;
			if ($save_files =~ /y/gi) {				
				my $output_filename = join ("_", @arguments);
				$output_filename = "Output/".$output_filename."_$HG_counter\_HGs";
	
				# 1) Save the groups of homology that match the query, with spp names
				#print "\nShowing first sequence names for first groups of homology:\n";
				#foreach (@good_homology_groups_spp_names) {
				#	print substr($_, 0, 160), "\n\n";
				#}
				unless (open (OUTPUT1, ">$output_filename"."_MCL_genes_IDs.out")) {
					print "Can't open file to save";
				}
				print OUTPUT1 @good_homology_groups_spp_names;
				close OUTPUT1;
				
				# 2) Save the columns from MCL parsed file for groups of homology that match the query
				#print "\nShowing columns for the few first groups of homology:\n";
				#print "Species: Cfra Cowc Saro Aque Emue Mlei Tadh Aequ Nvec Hvul Hmia Bcal Aric Aste Dcar Rsoc Rsor Rota Mlig Pcro Smed Lana Paus Bner Pesc Dgyr Lluy Pech Ctel Eand Mvul Hrob Hman Hmed Llon Hhan Obim Spha Dpol Mmer Mcal Myes Cgig Lgig Ppel Hrub Gmag Gaeg Mcor Acal Echl Bgla Aful Cuni Oida Pcau Tuco Psam Hali Ppac Dcor Cele Rhab Anan Pred Pdav Hexe Rvar Ebro Nstr Lpol Popi Dpte Dsil Cscu Ptep Gmae Hhol Smar Rimm Dste Eaff Tcal Ppol Pchi Hame Pcla Cqua Anas Hazt Dmag Dgal Dpul Dapu Tqin Caug Iele Znev Tcas Dmel Skow Anja Ajap Spur Arub Blan Cint Ebur Cmil Leri Cpla Drer Sorb Bpec Pmag Pmod Onil Olat Trub Lcha Pann Rbiv Gser Muni Pwal Bbom Xlae Xtro Sbom Bbuf Npar Rtem Emac Acar Ctig Psin Deco Care Cmyd Cser Pmeg Cpic Tcar Amis Asin Cpor Ggan Ggal Tgut Oana Mdom Hsap Clup Ttru\n\n";
				#foreach (@good_homology_groups_columns_parsed) {
				#	print $_;
				#}
				unless (open (OUTPUT2, ">$output_filename"."_MCL_columns_parsed.out")) {
					print "Can't open file to save";
				}
				print OUTPUT2 "HG\tCfra\tCowc\tSaro\tAque\tEmue\tMlei\tTadh\tAequ\tNvec\tHvul\tHmia\tBcal\tAric\tAste\tDcar\tRsoc\tRsor\tRota\tMlig\tPcro\tSmed\tLana\tPaus\tBner\tPesc\tDgyr\tLluy\tPech\tCtel\tEand\tMvul\tHrob\tHman\tHmed\tLlon\tHhan\tObim\tSpha\tDpol\tMmer\tMcal\tMyes\tCgig\tLgig\tPpel\tHrub\tGmag\tGaeg\tMcor\tAcal\tEchl\tBgla\tAful\tCuni\tOida\tPcau\tTuco\tPsam\tHali\tPpac\tDcor\tCele\tRhab\tAnan\tPred\tPdav\tHexe\tRvar\tEbro\tNstr\tLpol\tPopi\tDpte\tDsil\tCscu\tPtep\tGmae\tHhol\tSmar\tRimm\tDste\tEaff\tTcal\tPpol\tPchi\tHame\tPcla\tCqua\tAnas\tHazt\tDmag\tDgal\tDpul\tDapu\tTqin\tCaug\tIele\tZnev\tTcas\tDmel\tSkow\tAnja\tAjap\tSpur\tArub\tBlan\tCint\tEbur\tCmil\tLeri\tCpla\tDrer\tSorb\tBpec\tPmag\tPmod\tOnil\tOlat\tTrub\tLcha\tPann\tRbiv\tGser\tMuni\tPwal\tBbom\tXlae\tXtro\tSbom\tBbuf\tNpar\tRtem\tEmac\tAcar\tCtig\tPsin\tDeco\tCare\tCmyd\tCser\tPmeg\tCpic\tTcar\tAmis\tAsin\tCpor\tGgan\tGgal\tTgut\tOana\tMdom\tHsap\tClup\tTtru\n";
				print OUTPUT2 @good_homology_groups_columns_parsed;
				close OUTPUT2;
				
				# 3) Now save in return format, one taxa per line
				#@good_homology_groups_spp_names = @good_homology_groups_spp_names_copy;
				my @gene_names_return = ();
				print "\nShowing first sequence names for groups of homology:\n";
				foreach (@good_homology_groups_spp_names) {
					chomp $_;
					my @homology_group = split (" ", $_);
					my $gene_names = parse_gene_names_return_format(@homology_group);
					push (@gene_names_return, $gene_names);
					print "\n",substr($gene_names, 0, 480),"...\n";
				}
				unless (open (OUTPUT3, ">$output_filename"."_MCL_annotated_genes.out")) {
					print "Can't open file to save";
				}
				print OUTPUT3 @gene_names_return;
				close OUTPUT3;
				
				# 4) Save the names of the taxa present in the groups of homology 
				my @taxa_names = ();
				foreach (@good_homology_groups_spp_names) {
					chomp $_;
					my @homology_group = split (" ", $_);
					push (@taxa_names, parse_taxa_names(@homology_group));
				}
				#print "\nShowing the labels of all the taxa for each homology group:\n";
				#foreach (@taxa_names) {
				#	print "$_";
				#	}
				unless (open (OUTPUT4, ">$output_filename"."_taxa_names.out")) {
					print "Can't open file to save";
				}
				print OUTPUT4 @taxa_names;
				close OUTPUT4;
				print "\nResults saved to files $output_filename\n";
			}
		}
	}
} until ($user_input =~ /^\s*$|^exit|^quit/i);

#Untie hash files
untie %MCL_out;
untie %MCL_columns_parsed;

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
    unless (-e $filename) {print "File $filename does NOT exists\n" and exit;}
    unless (-f $filename) {print "File $filename is NOT a flat file\n" and exit;}
    unless (-s $filename) {print "File $filename is empty\n" and exit;}
    unless (open (FH, $filename)) {print "File $filename can not be opened\n" and exit;}
    close FH;
    return;
}

#MCL_COLUMNS_PARSED_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_columns_parsed_to_hash {
	my ($filename) = @_;
	my %hash;
	my @line = '';
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {	#Parse file one line at a time
		$line_counter++;
		chomp $line;
		$line =~ s/^\d*\s//;	#To remove first digits, which indicate the group of homology number/ID
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#MCL_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_output_to_hash {
	my ($filename) = @_;
	my %hash;
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {   			#Parse file one line at a time
		$line_counter++;
		chomp $line;
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#PARSE_GENE_NAMES_TEXT_FORMAT: subrout to extract the gene names from the annotated genomes for the list of groups of homology
# producing an output in text format
sub parse_gene_names_text_format {
	my (@homology_group) = @_;
	my @list_annotated_genomes = qw /HG_ Cowc Saro Aque Tadh Nvec Hvul Bcal Aste Dcar Rsoc Rsor Rota Mlig Lana Bner Dgyr Ctel Hrob Obim Spha Dpol Mmer Mcal Myes Cgig Lgig Hrub Gaeg Acal Echl Bgla Cuni Pcau Hali Ppac Cele Hexe Nstr Lpol Dpte Dsil Cscu Ptep Dste Eaff Tcal Ppol Pchi Hame Pcla Cqua Anas Hazt Dmag Dgal Dpul Dapu Tqin Iele Znev Tcas Dmel Skow Anja Ajap Spur Arub Cint Ebur Cmil Leri Cpla Drer Sorb Bpec Pmag Onil Olat Trub Lcha Pann Rbiv Gser Muni Bbom Xlae Xtro Sbom Bbuf Npar Rtem Emac Acar Ctig Psin Deco Care Cmyd Cser Pmeg Cpic Tcar Amis Asin Cpor Ggan Ggal Tgut Oana Mdom Hsap Clup Ttru/;
	#@list_annotated_genomes = reverse @list_annotated_genomes;
	my $results;

	foreach my $annotated_taxa (@list_annotated_genomes) {
		#print "\$gene: $gene\n";
		foreach my $gene (@homology_group) {
			#print "\$annotated_taxa: $annotated_taxa ";
			if ($gene =~ /$annotated_taxa/i) {
				#print color ("red"), "Match!\n", color ("reset");
				$results .= "$gene ";
			}
		}
	}
	#print "\@results in sub: $results\n";
	if ($results =~ /^HG_\d*\s*$/) {
		$results .= "This group of homology does not contain any annotated genome";
	}
	#print "Results @results\n";
	#@results = sort @results;
	$results .= "\n";
	return $results;
}

#PARSE_GENE_NAMES_RETURN_FORMAT: subrout to extract the gene names from the annotated genomes for the list of groups of homology
# producing an output in a format with returns
sub parse_gene_names_return_format {
	my (@homology_group) = @_;
	my @list_annotated_genomes = qw /HG_ Cowc Saro Aque Tadh Nvec Hvul Bcal Aste Dcar Rsoc Rsor Rota Mlig Lana Bner Dgyr Ctel Hrob Obim Spha Dpol Mmer Mcal Myes Cgig Lgig Hrub Gaeg Acal Echl Bgla Cuni Pcau Hali Ppac Cele Hexe Nstr Lpol Dpte Dsil Cscu Ptep Dste Eaff Tcal Ppol Pchi Hame Pcla Cqua Anas Hazt Dmag Dgal Dpul Dapu Tqin Iele Znev Tcas Dmel Skow Anja Ajap Spur Arub Cint Ebur Cmil Leri Cpla Drer Sorb Bpec Pmag Onil Olat Trub Lcha Pann Rbiv Gser Muni Bbom Xlae Xtro Sbom Bbuf Npar Rtem Emac Acar Ctig Psin Deco Care Cmyd Cser Pmeg Cpic Tcar Amis Asin Cpor Ggan Ggal Tgut Oana Mdom Hsap Clup Ttru/;
	#@list_annotated_genomes = reverse @list_annotated_genomes;
	my $results;
	my $flag = 0;

	foreach my $annotated_taxa (@list_annotated_genomes) {
		$flag = 0;
		#print "\$gene: $gene\n";
		foreach my $gene (@homology_group) {
			#print "\$annotated_taxa: $annotated_taxa ";
			if ($gene =~ /$annotated_taxa/i) {
				#print color ("red"), "Match!\n", color ("reset");
				$results .= "$gene\t";
				$flag = 1;
			}
		}
		if ($flag == 1) { $results .= "\n"; }
	}
	
	#print "\@results in sub: $results\n";
	if ($results =~ /^HG_\d*\s*$/) {
		$results .= "This group of homology does not contain any annotated genome\n";
	}
	#print "Results @results\n";
	#@results = sort @results;
	#$results .= "\n";
	return $results;
}

#PARSE_TAXA_NAMES: subrout to extract the taxa names from the list of groups of homology
sub parse_taxa_names {
	my (@homology_group) = @_;
	my %taxons;
	my $results;
	my $group_ID = splice (@homology_group, 0 ,1);

	foreach my $taxon (@homology_group) {
		$taxon =~ /^.{4}/;
		$taxons{$&} = '';
	}
	my @keys = sort keys %taxons;
	$results = join ("\t", @keys);
	$results = "$group_ID "."\t$results"."\n";
	return $results;
}
#PRINT_EXAMPLES: subroutine to print some examples to user
sub print_examples {
	print '
Search the MCL groups of homology by gene content.
The search is based in the gene names provided by the annotated genomes, such as Homo sapiens.
This is a very basic text search, and searches independently for each of the terms introduce separated by spaces
Search is case insensitive.

Examples:
	Hox1
	wingless
	TGF-Beta
	
More than one gene name can be used in a single serach, separated by spaces:
	Hox1 Cdx Hox10
	wingless wnt2 wnt4 wnt4
	TGF_Beta BMP activin
';
	return;
}