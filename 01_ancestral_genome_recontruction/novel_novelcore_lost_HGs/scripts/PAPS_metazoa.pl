#!/usr/bin/perl -w
#Jordi Paps, Oxford 2015; with help of Patrick Gemmell (Oxford) for some subroutines
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
print "Species: Cfra Cowc Saro Aque Emue Mlei Tadh Aequ Nvec Hvul Hmia Bcal Aric Aste Dcar Rsoc Rsor Rota Mlig Pcro Smed Lana Paus Bner Pesc Dgyr Lluy Pech Ctel Eand Mvul Hrob Hman Hmed Llon Hhan Obim Spha Dpol Mmer Mcal Myes Cgig Lgig Ppel Hrub Gmag Gaeg Mcor Acal Echl Bgla Aful Cuni Oida Pcau Tuco Psam Hali Ppac Dcor Cele Rhab Anan Pred Pdav Hexe Rvar Ebro Nstr Lpol Popi Dpte Dsil Cscu Ptep Gmae Hhol Smar Rimm Dste Eaff Tcal Ppol Pchi Hame Pcla Cqua Anas Hazt Dmag Dgal Dpul Dapu Tqin Caug Iele Znev Tcas Dmel Skow Anja Ajap Spur Arub Blan Cint Ebur Cmil Leri Cpla Drer Sorb Bpec Pmag Pmod Onil Olat Trub Lcha Pann Rbiv Gser Muni Pwal Bbom Xlae Xtro Sbom Bbuf Npar Rtem Emac Acar Ctig Psin Deco Care Cmyd Cser Pmeg Cpic Tcar Amis Asin Cpor Ggan Ggal Tgut Oana Mdom Hsap Clup Ttru\n\n";
my @MCL_columns_parsed_keys = keys %MCL_columns_parsed;
foreach (@MCL_columns_parsed_keys[0...4]){
	print "$_ => $MCL_columns_parsed{$_}\n";
} print "...\n";

#Create hash from the clade definition in subrout CLADES
my %spp;
my @value_list = ();
%spp = hash_spp();

#Calculate how many columns (species/terminal tip) are in the hash
walk_hash (\%spp, \@value_list);
my @total_columns = scalar @value_list;
print "\nNumber of species: ", (scalar @value_list), "\n";

#Ask user for clade/spp to check genes taxonomic distribution, perform search
my $user_input = '';	#To store the searching criteria fom user
my @arguments = ();	#To store the split of the searching criteria
my @search= ();  	#To store the taxa and options of the searching criteria
my $taxa  = '';		#To store the taxa from @search
my $option = '';	#To store options from @search (present, absent, minus, atleast,...)
my @columns = ();	#To store columns to search
my $final_search;	#To store columns to search plus the options (present, absent, minus, atleast,...)
my %final_searches = ();#To store ALL searches columns and options (present, absent, minus, atleast,...)
my @true_flags = (); 	#To store the flags that indicate if a homology group/line passes all the queries checks
my @outgroup = ();	#To store the columns left over at the end after extracting all the ingroups columns
my @good_homology_groups_spp_names = ();	#To store the groups of homology fullfilling the search criteria
my @good_homology_groups_spp_names_copy = ();	#Backup
my @good_homology_groups_columns_parsed = ();	#To store the columns of the groups of homology fullfilling the search criteria

OUTER: do {
	print "\nPlease, enter name of the clade/species to search (\"example\" for some samples, \"tree\" to print the evolutionary tree, Enter to exit):\n";
	$user_input = <STDIN>;
	chomp $user_input;
	unless ($user_input =~ /^\s*$|^exit|^quit/i) {			#Unless user wants to exit...
		if ($user_input =~ /example|help/i) {
			print_examples();			#Print examples of commands if user requests it
		} elsif ($user_input =~ /tree/i) {
			print_hash_colors(%spp);		#Print examples of commands if user requests it
		} else {
			#Here the real search starts, emptying variables for next loop and parsing the user input
			%final_searches = ();			#Empty the hash containing all the search conditions
			@good_homology_groups_spp_names = ();	#Empty the hash containing the results from previous search
			@good_homology_groups_columns_parsed = (); #Empty the hash containing the columns of the groups of homology fullfilling the search criteria
			@arguments = split (" ", $user_input);	#Decompose the user input in different arguments entered by user, each taxa query delimited by a space (" ")
			#print "\@arguments: @arguments\n";			
			foreach (@arguments) {
				#print "\$_ in \@arguments: $_\n";
				@search = split ("-", $_);	#Decompose each argument into taxa (first item) and options (second item)
				#print "\@search: @search\n";
				$taxa = $search[0];
				$option = $search[1];
				unless (defined $option && $option =~ /pre|all|abs|none|min|but|atl|onl|jus/gi){
						print "Taxa $taxa is missing valid options (present, absent, minus, atleast, only,...)\n";
						goto OUTER;
				}
				if ($taxa =~ /^out|^rest|^other/i) {		#We store outgroup conditions in the hash, as "outgroup" does not exists in the species hash of hashes
					$final_searches {$taxa."_".$option} = "outgroup_".$option;       #Store ALL the searches in a hash, keys are taxa and values the columns and options
				} else {
					@columns = obtain_taxa_columns (\$taxa, \%spp, \%MCL_columns_parsed); #Obtain the columns belonging to each taxa
					#print "Columns that will be inspected: @columns\n";
					if (scalar @columns == 0 ){
						print "Taxa $taxa not found in taxa list\n";
						goto OUTER;
					}
					unless (defined $search[1] && $search[1] =~ /pre|all|abs|none|min|but|atl|onl|jus/gi){
						print "Taxa $taxa is missing valid options (present, absent, minus, atleast, only)\n";
						goto OUTER;
					}
					$final_search = join ("_", @columns, $option); 		#Join the columns and the options to store them later in %final_search
					$final_searches {$taxa."_".$option} = $final_search;	#Store ALL the searches in a hash, keys are taxa and values the columns and options
				}	
			}
			
			#Read the hash containing the file with the MCL output parsed to columns, check the search criteria
			my @keys_MCL_columns_parsed = keys %MCL_columns_parsed;			
			#print "Colums to be inspected: @keys_MCL_columns_parsed\n";
			my $HG_counter = 0;
			foreach my $homology_group (@keys_MCL_columns_parsed){ 				#For each line of the MCL parsed file...
				#print "\n", '$homology_group: ', $homology_group, "\n";
				@true_flags = ();							#Empty the flags that indicate if a homology group/line passes all the queries checks
				my @MCL_columns = split (" ", $MCL_columns_parsed{$homology_group});	#Explode the string of the homology group to columns				
				@outgroup = @MCL_columns;						#Storing all columns in @outgroups, later the ingroups columns will be spliced out o this array
				my @keys = keys %final_searches;
				#print 'keys in %final_searches: ', "@keys" , "\n";
				#print "MCL columns_parsed:\n@MCL_columns\n";
				foreach my $query (keys %final_searches) {				#Now go query from query stored in %final_searches
					#print '$query ', $query,"\n";
					unless ($query =~ /outg/i) {					#Leave outgroup check for the end
						my @query = split ("_", $final_searches{$query});	#Explode the string of the columns to check
						#print "\@query: @query ";
						my $condition = pop @query;				#Save the options (present, absent, minus, atleast,...)
						#print "\$condition: $condition \n";
						@query = @query;				
						#print "\@query columns:  @query\n";
						my $total_taxa_in_clade = scalar @query;		#Stores total number of taxa in queried clade, used to later check options
						my $sum = 0;
						#print "Columns contens: ";
						foreach my $query_column (@query) {			#For each column corresponding to the taxa in the query, we extract the value of the column in MCL file
							#print "$MCL_columns[$query_column] ";							
							if ($MCL_columns[$query_column] != 0) {
								$sum++;
							}
							splice (@outgroup, $query_column, 1, "X"); 	#Remove column from outgroup array, to check later the remaining columns
						}
						#print "\nTotal taxa in clade:$total_taxa_in_clade\n";
						#print "Taxa of clade in this MCL: $sum\n";
						#Check if the group of homology matches the conditions requested by user for this clade
						my $check = check_conditions($condition, $total_taxa_in_clade, $sum);
						if ($check eq "wrong") {
							goto OUTER;
						} else {
							push (@true_flags, $check);
						}
					}
				}
				#OUTGROUP CHECK
				foreach (keys %final_searches) {				#Now we check if user asked for outgroup conditions and is stored in %final_searches
					if ($_ =~ /outg/i) {
						#print "OUTGROUP\n";
						#print "\@outgroup:\n@outgroup\n";
						#print "\@MCL_columns:\n@MCL_columns\n";
						my @query = split ("_", $final_searches{$_});	#Explode the string of the columns to check
						my $condition = pop @query;				#Save the options (present, absent, minus, atleast,...)
						#print '$condition in outgroup: ', $condition,"\n";
						@outgroup = grep { $_ ne "X" } @outgroup;
						#print 'MCL contens in outgroup: ', "\n@outgroup","\n";
						my $total_taxa_in_clade = scalar @outgroup;		#Stores total number of taxa in queried clade, used to later check options
						#print 'Total taxa in outgroup: ', $total_taxa_in_clade, "\n";
						my $sum = 0;
						#print "\$MCL_columns[\$query_column] in outgroup:\n";
						foreach (@outgroup) {			#For each column corresponding to the taxa in the query, we extract the value of the column in MCL file
							#print "$_ ";
							if ($_ != 0) {
								$sum++;
							}
						}
						#print 'Outgroup in this MCL: ', $sum, "\n";
						my $check = check_conditions($condition, $total_taxa_in_clade, $sum);
							if ($check eq "wrong") {
								goto OUTER;
							} else {
								push (@true_flags, $check);
							}
					}
				}
				#print "\@true_flags: @true_flags\n";
				#Check if all flags are true, then store the group of homology in %results
				my $flags = 0;
				
				foreach (@true_flags) {
					if ($_ !~ 'true' ) {
						$flags++;
					}					
				}
				if ($flags == 0) {
					push (@good_homology_groups_spp_names, "$homology_group\t$MCL_out{$homology_group}\n");
					push (@good_homology_groups_columns_parsed, "$homology_group\t$MCL_columns_parsed{$homology_group}\n");
					$HG_counter++;
				}				
			}
			@good_homology_groups_spp_names_copy = @good_homology_groups_spp_names;
			print "\nNumber of groups of homology found: $HG_counter \n";
			if ($HG_counter == 0) { goto OUTER; };
			
			#Save 4 different output files
			print "\nDo you want to see results (gene names, MCL groups, etc) and save them in files? Yes/No\n";
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
	foreach my $line (@file) {   	#Parse file one line at a time
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

#WALK_HASH: subroutine to traverse the hash of hashes, modified after in http://www.perlmonks.org/?node_id=116162
sub walk_hash { 
my ($hash, $value_list, $key_list) = @_;
while (my ($key, $value) = each %$hash) {
	push @$key_list, $key;
	if (ref($value) eq 'HASH') {
		walk_hash ($value, $value_list, $key_list);
	} else {
		push @$value_list, $value;
	}
	pop @$key_list;
	}
}

#QUERY: intermediate subroutine that sends the taxa search to the to recoursive subroutines PG_WALK and PRINT_EVERYTHING
sub obtain_taxa_columns { 
	my ($taxa, $spp, $MCL_parsed) = @_; 
	my @results = ();
	
	#Send query and hash      

	pg_walk (\%$spp, [], $$taxa, \@results);
	@results = sort {$a <=> $b} @results;
	return @results; 
}

#PG_WALK: subroutine to traverse the hash of hashes till finding the queried label
#Hat tip to Patrick Gemmell (Oxford) for the help, he modified the subroutine
#found in http://www.perlmonks.org/?node_id=116162
sub pg_walk {
	my ($hash, $key_list, $query_text, $localresults) = @_;
	while (my ($key, $value ) = each %$hash) {
		push @$key_list, $key;
		if ($key =~ /^$query_text/gi) {
			print "Taxa that will be searched: $key\n";
			print_everything($value , $localresults);
		} else {
			if (ref($value ) eq 'HASH') {
				pg_walk($value , $key_list, $query_text, $localresults);
			}
		}
		pop @$key_list;
    }
}

#PRINT_EVERYTHING: subroutine to traverse the hash of hashes to print the values corresponding to the query in sub PG_WALK (see above).
#Hat tip to Patrick Gemmell (Oxford) for the help, he modified the subroutine
#found in http://www.perlmonks.org/?node_id=116162
sub print_everything {
	my ($hash, $localresults, $key_list) = @_;
	while (my ($key, $value) = each %$hash) {
		push @$key_list, $key;
		if (ref($value) eq 'HASH') {
			print_everything($value, $localresults, $key_list);
		} else {
			push @$localresults, $value;
		}
		pop @$key_list;
	}
}

#CHECK_CONDITIONS: subroutine to check the conditions of presence/absence specified by the user, returns TRUE or FALSE
sub check_conditions {
	my ($condition, $total_taxa_in_clade, $sum) = @_;
	#Perform the check according to the different options (present, absent, minus, atleast, only...
	if ($condition =~ /pre|all/gi) {
		if ($sum == $total_taxa_in_clade) {
			#print "True all\n";
			return "true";
		} else {
			#print "False all\n";
			return "false";
		}
	} elsif ($condition =~ /abs|none/gi) {
		if ($sum == 0) {
			#print "True none\n";
			return "true";
		} else {
			#print "False none\n";
			return "false";
		}
	} elsif ($condition =~ /atl\D*/gi ) {
		$condition =~ s/$&//g;
		#print "\$condition: $condition\n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"at least\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum >= $condition) {
			#print "True atleast $condition\n";
			return "true";
		} else {
			#print "Atleast $condition false\n";
			return "false";
		}
	} elsif ($condition =~ /min\D*|but\D*/gi ) {
		$condition =~ s/$&//g;							
		#print "\$condition: $condition \n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"minus\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum == ($total_taxa_in_clade - $condition)) {
			#print "True minus $condition\n";
			return "true";
		} else {
			#print "Minus $condition false\n";
			return "false";
		}
	} elsif ($condition =~ /onl\D*|jus\D*/gi ) {
		$condition =~ s/$&//g;								
		#print "\$condition: $condition \n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"only\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum == $condition) {
			#print "True only $condition\n";
			return "true";
		} else {
			#print "Only $condition false\n";
			return "false";
		}
	}
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
Clade/species names can be truncated, but the start of the clade name should match the table printed above.
Search is case insensitive.

Some search examples (first 4 digits in examples stand for rest of taxa, the other 4 for ingroup):	
  "Vertebrata-present" => genes found in ALL vertebrate species, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1111
  "Vertebrata-present Rest-absent" => genes found in ALL vertebrate species, absent in other clades/rest of taxa
	Rest of taxa 0000 Ingroup 1111
  "Vertebrata-present Rest-present" => genes found in ALL vertebrate species, present in other clades/rest of taxa
	Rest of taxa 1111 Ingroup 1111
  "Vertebrata-absent Rest-present" => genes found in rest of taxa species, absent in Vertebrata
	Rest of taxa 1111 Ingroup 0000
  "Homo-present Mus-present Rest-absent" => genes only found in humans and mice. Species can be specified one by one.

The number of species presenting/missing for a gene can be fine-tuned with minus#, atleast#, only# for both ingroup and rest of taxa:
  "Vertebrata-minus1" => found in ALL vertebrate species but one, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1110 / 1101 / 1011 / 0111
  "Vertebrata-minus2 Rest-minus1" => genes found in ALL vertebrate species but one, absent in other clades/rest of taxa
	Rest of taxa 1110 / 1101 / 1011 / 0111 Ingroup 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Vertebrata-atleast1 Rest-atleast1" => genes found in at least 1 vertebrate species and 1 rest of taxa species
	Rest of taxa 1000 / 1100 / 1110 / 1111 / 1010 / 1011 / 1001 / 1101 / 0110 / 0111 etc.
	Ingroup  1000 / 1100 / 1110 / 1111 / 1010 / 1011 / 1001 / 1101 / 0110 / 0111 etc.
  "Vertebrata-only3" => return genes found in just 3 vertebrate species, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1110 / 1101 / 1011 / 0111
  
Different criteria can be combined in a single search:
  "Vertebrata-minus1 Echinodermata-atleast2" => genes found in ALL vertebrate species but one, AND present in at least two echinoderms, absent/present in other clades/rest of taxa
	Rest of taxa ???? Vertebrata 1110 / 1101 / 1011 / 0111 Echinodermata 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Vertebrata-atleast2 Urochordata-atleast2" => genes found in 2 or more vertebrate species OR 2 or more urochordates, independently if they are present/absent in other clades/rest of taxa
	Rest of taxa ???? Vertebrata 1100 / 1010 / 1001 / 0110 / 0101 / 0011 Urochordata 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Nematoda-absent Platyhelminthes-absent Rest-present" => genes found in clades/rest of taxa, absent (convergently lost) in round worms and flatworms
	Rest of taxa 1111 Nematoda 0000 Platyhelminthes 0000
	
Carefull with nested taxa!!! Start with the greater group taking into account the conditions for the smaller group:
  To find genes in ALL chordates but missing only in humans => "Chordata-minus1 Hsap-absent"
  To find genes in ALL chordates but missing only in vertebrates => "Chordata-minus5 Vertebrata-absent"
  To find genes in at least one clade of chordates, but missing only in vertebrates => "Cephalocordata-atleast1 Urochordata-atleast1 Vertebrata-absent"
  ';
	return;
}

#HASH_SPP: a subroutine to define define the hash of hashes containing of all the clades and spp included, put here to not clutter the main program.
#Each species is assigned a numeric value, same as the column they occupy in the parsed MCL output, so Crei is the first column and Hsap occupies the last column.
#Thus, when user asks for a group, these values can be used as index to lookup lines/arrays.
#Each element shoud have the same number of levels, or the subrout "PRINT_HASH_COLORS" won't work.
sub hash_spp {
	my (%spp) = ();

#All spp, one by one
###   domain # node_a # node_b # node_c # phylum_a # phylum_b # phylum_c # phylum_d # phylum_e # phylum_f # phylum_g # phylum_h # phylum_i # phylum_j # phylum_k # subphylum_a # subphylum_b # subphylum_c # subphylum_d # superclass # class_a # class_b # class_c # class_d # class_e # class_f # class_g # subclass_a # subclass_b # subclass_c # order_a # order_b # order_c # order_d # suborder # infraorder # family_a # family_b # genus # species
$spp {'Eukaryota'}{'Ichthyosporea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Creolimax_fragrantissima_(Cfra)'}{''} = 0;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Filasterea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Capsaspora_owczarzaki_(Cowc)'}{''} = 1;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Choanoflagellatea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Salpingoeca_rosetta_(Saro)'}{''} = 2;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Porifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Amphimedon'}{'Amphimedon_queenslandica_(Aque)'}{''} = 3;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Porifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Ephydatia'}{'Ephydatia_muelleri_(Emue)'}{''} = 4;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Ctenophora'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Mnemiopsis'}{'Mnemiopsis_leidyi_(Mlei)'}{''} = 5;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Placozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Trichoplax'}{'Trichoplax_adhaerens_(Tadh)'}{''} = 6;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Cnidaria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Anthozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Actinia'}{'Actinia_equina_(Aequ)'}{''} = 7;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Cnidaria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Anthozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Nematostella'}{'Nematostella_vectensis_(Nvec)'}{''} = 8;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Cnidaria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Hydrozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Hydra'}{'Hydra_vulgaris_(Hvul)'}{''} = 9;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Xenacoelomorpha'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Hofstenia'}{'Hofstenia_miamia_(Hmia)'}{''} = 10;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Monogononta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Ploima'}{''}{''}{''}{''}{''}{''}{''}{'Brachionus'}{'Brachionus_calyciflorus_(Bcal)'}{''} = 11;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloidea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Adinetida'}{''}{''}{''}{''}{''}{''}{''}{'Adineta_r'}{'Adineta_ricciae_(Aric)'}{''} = 12;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloidea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Adinetida'}{''}{''}{''}{''}{''}{''}{''}{'Adineta_s'}{'Adineta_steineri_(Aste)'}{''} = 13;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloidea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloida'}{''}{''}{''}{''}{''}{''}{''}{'Didymodactylos'}{'Didymodactylos_carnosus_(Dcar)'}{''} = 14;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloidea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloida'}{'Bdelloida_a'}{''}{''}{''}{''}{''}{''}{'Rotaria_c'}{'Rotaria_socialis_(Rsoc)'}{''} = 15;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloidea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloida'}{'Bdelloida_a'}{''}{''}{''}{''}{''}{''}{'Rotaria_r'}{'Rotaria_sordida_(Rsor)'}{''} = 16;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Rotifera'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloidea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bdelloida'}{'Bdelloida_a'}{''}{''}{''}{''}{''}{''}{'Rotaria_p'}{'Rotaria_sp._Silwood1_(Rota)'}{''} = 17;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Platyhelminthes'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Macrostomum'}{'Macrostomum_lignano_(Mlig)'}{''} = 18;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Platyhelminthes'}{'Platyhelminthes_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Prostheceraeus'}{'Prostheceraeus_crozieri_(Pcro)'}{''} = 19;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Platyhelminthes'}{'Platyhelminthes_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Schmidtea'}{'Schmidtea_mediterranea_(Smed)'}{''} = 20;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophophorata'}{'Lophophorata_a'}{'Brachiopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Lingula'}{'Lingula_anatina_(Lana)'}{''} = 21;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophophorata'}{'Lophophorata_a'}{'Phoronida'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Phoronis'}{'Phoronis_australis_(Paus)'}{''} = 22;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophophorata'}{''}{'Bryozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bugula'}{'Bugula_neritina_(Bner)'}{''} = 23;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Sipuncula'}{''}{''}{''}{''}{''}{''}{''}{'Phascolosoma'}{'Phascolosoma_esculenta_(Pesc)'}{''} = 24;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Eunicida'}{''}{''}{''}{''}{''}{''}{''}{'Dimorphilus'}{'Dimorphilus_gyrociliatus_(Dgyr)'}{''} = 25;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Sabellida'}{''}{''}{''}{''}{''}{''}{''}{'Lamellibrachia'}{'Lamellibrachia_luymesi_(Lluy)'}{''} = 26;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Sabellida'}{''}{''}{''}{''}{''}{''}{''}{'Paraescarpia'}{'Paraescarpia_echinospica_(Pech)'}{''} = 27;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{'Annelida_c'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Capitellidae'}{''}{'Capitella'}{'Capitella_teleta_(Ctel)'}{''} = 28;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{'Annelida_c'}{''}{''}{'Clitellata'}{''}{''}{''}{''}{''}{''}{'Oligochaeta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Eisenia'}{'Eisenia_andrei_(Eand)'}{''} = 29;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{'Annelida_c'}{''}{''}{'Clitellata'}{''}{''}{''}{''}{''}{''}{'Oligochaeta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Metaphire'}{'Metaphire_vulgaris_(Mvul)'}{''} = 30;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{'Annelida_c'}{''}{''}{'Clitellata'}{''}{''}{''}{''}{''}{''}{'Hirudinea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Helobdella'}{'Helobdella_robusta_(Hrob)'}{''} = 31;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{'Annelida_c'}{''}{''}{'Clitellata'}{''}{''}{''}{''}{''}{''}{'Hirudinea'}{''}{''}{'Arhynchobdellida'}{''}{''}{''}{''}{''}{''}{''}{'Hirudinaria'}{'Hirudinaria_manillensis_(Hman)'}{''} = 32;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Annelida'}{'Annelida_a'}{'Annelida_b'}{'Annelida_c'}{''}{''}{'Clitellata'}{''}{''}{''}{''}{''}{''}{'Hirudinea'}{''}{''}{'Arhynchobdellida'}{''}{''}{''}{''}{''}{''}{''}{'Hirudo'}{'Hirudo_medicinalis_(Hmed)'}{''} = 33;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Lophotrochozoa_d'}{'Nemertea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Lineus'}{'Lineus_longissimus_(Llon)'}{''} = 34;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Polyplacophora'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Hanleya'}{'Hanleya_hanleyi_(Hhan)'}{''} = 35;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Cephalopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Octopus'}{'Octopus_bimaculoides_(Obim)'}{''} = 36;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Cephalopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Sepia'}{'Sepia_pharaonis_(Spha)'}{''} = 37;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Bivalvia'}{'Bivalvia_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Dreissena'}{'Dreissena_polymorpha_(Dpol)'}{''} = 38;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Bivalvia'}{'Bivalvia_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Mercenaria'}{'Mercenaria_mercenaria_(Mmer)'}{''} = 39;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Bivalvia'}{'Bivalvia_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Mytilus'}{'Mytilus_californianus_(Mcal)'}{''} = 40;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Bivalvia'}{'Bivalvia_b'}{'Bivalvia_c'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Mizuhopecten'}{'Mizuhopecten_yessoensis_(Myes)'}{''} = 41;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Bivalvia'}{'Bivalvia_b'}{'Bivalvia_c'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Crassostrea'}{'Crassostrea_gigas_(Cgig)'}{''} = 42;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Patellogastropoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Lottia'}{'Lottia_gigantea_(Lgig)'}{''} = 43;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Patellogastropoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Patella'}{'Patella_pellucida_(Ppel)'}{''} = 44;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Vetigastropoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Haliotis'}{'Haliotis_rubra_(Hrub)'}{''} = 45;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Vetigastropoda'}{''}{''}{''}{''}{''}{''}{''}{'Trochidae'}{''}{'Gibbula'}{'Gibbula_magus_(Gmag)'}{''} = 46;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Vetigastropoda'}{''}{''}{''}{''}{''}{''}{''}{'Trochidae'}{''}{'Gigantopelta'}{'Gigantopelta_aegis_(Gaeg)'}{''} = 47;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Caenogastropoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Marisa'}{'Marisa_cornuarietis_(Mcor)'}{''} = 48;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Heterobranchia'}{''}{''}{''}{''}{''}{''}{''}{''}{'Aplysia'}{'Aplysia_californica_(Acal)'}{''} = 49;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Heterobranchia'}{'Heterobranchia_a'}{''}{''}{''}{''}{''}{''}{''}{'Elysia'}{'Elysia_chlorotica_(Echl)'}{''} = 50;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Heterobranchia'}{'Heterobranchia_a'}{'Heterobranchia_b'}{''}{''}{''}{''}{''}{''}{'Biomphalaria'}{'Biomphalaria_glabrata_(Bgla)'}{''} = 51;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Heterobranchia'}{'Heterobranchia_a'}{'Heterobranchia_b'}{'Stylommatophora'}{''}{''}{''}{''}{''}{'Achatina'}{'Achatina_fulica_(Aful)'}{''} = 52;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Heterobranchia'}{'Heterobranchia_a'}{'Heterobranchia_b'}{'Stylommatophora'}{'Stylommatophora_a'}{''}{''}{''}{''}{'Candidula'}{'Candidula_unifasciata_(Cuni)'}{''} = 53;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Lophotrochozoa'}{'Lophotrochozoa_a'}{'Lophotrochozoa_b'}{'Lophotrochozoa_c'}{'Mollusca'}{''}{''}{''}{''}{''}{''}{'Mollusca_a'}{'Mollusca_b'}{'Gastropoda'}{''}{''}{''}{''}{'Gastropoda_a'}{'Gastropoda_b'}{'Heterobranchia'}{'Heterobranchia_a'}{'Heterobranchia_b'}{'Stylommatophora'}{'Stylommatophora_a'}{''}{''}{''}{''}{'Oreohelix'}{'Oreohelix_idahoensis_(Oida)'}{''} = 54;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Priapulida'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Priapulus'}{'Priapulus_caudatus_(Pcau)'}{''} = 55;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Priapulida'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Tubiluchus'}{'Tubiluchus_corallicola_(Tuco)'}{''} = 56;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadorea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Plectida'}{''}{''}{''}{''}{''}{''}{''}{'Plectus'}{'Plectus_sambesii_(Psam)'}{''} = 57;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadorea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Rhabditina'}{''}{''}{''}{'Halicephalobus'}{'Halicephalobus_sp._NKZ332_(Hali)'}{''} = 58;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadorea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Rhabditina'}{''}{'Rhabditina_a'}{''}{'Pristionchus'}{'Pristionchus_pacificus_(Ppac)'}{''} = 59;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadorea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Rhabditina'}{''}{'Rhabditina_a'}{'Rhabditidae'}{'Diploscapter'}{'Diploscapter_coronatus_(Dcor)'}{''} = 60;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadorea'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Rhabditina'}{''}{'Rhabditina_a'}{'Rhabditidae'}{'Caenorhabditis'}{'Caenorhabditis_elegans_(Cele)'}{''} = 61;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Rhabditophanes'}{'Rhabditophanes_sp._(Rhab)'}{''} = 62;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria_a'}{''}{'Acrobeloides'}{'Acrobeloides_nanus_(Anan)'}{''} = 63;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria_a'}{'Panagrolaimidae'}{'Panagrellus'}{'Panagrellus_redivivus_(Pred)'}{''} = 64;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Nematoda'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Chromadoria_a'}{'Panagrolaimidae'}{'Panagrolaimus'}{'Panagrolaimus_davidi_(Pdav)'}{''} = 65;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Tardigrada'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Hypsibius'}{'Hypsibius_exemplaris_(Hexe)'}{''} = 66;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Tardigrada'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Ramazzottius'}{'Ramazzottius_varieornatusa_(Rvar)'}{''} = 67;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Onychophora'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Epiperipatus'}{'Epiperipatus_broadwayi_(Ebro)'}{''} = 68;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Nymphon'}{'Nymphon_striatum_(Nstr)'}{''} = 69;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{'Euchelicerata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Xiphosura'}{''}{''}{''}{''}{''}{''}{''}{'Limulus'}{'Limulus_polyphemus_(Lpol)'}{''} = 70;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{'Euchelicerata'}{''}{''}{''}{'Arachnida'}{'Arachnida_a'}{''}{''}{''}{''}{''}{''}{''}{''}{'Opiliones'}{''}{''}{''}{''}{''}{''}{''}{'Phalangium'}{'Phalangium_opilio_(Popi)'}{''} = 71;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{'Euchelicerata'}{''}{''}{''}{'Arachnida'}{'Arachnida_a'}{''}{''}{''}{''}{''}{'Acari'}{''}{''}{'Acariformes'}{''}{''}{''}{''}{''}{''}{''}{'Dermatophagoides'}{'Dermatophagoides_pteronyssinus_(Dpte)'}{''} = 72;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{'Euchelicerata'}{''}{''}{''}{'Arachnida'}{'Arachnida_a'}{''}{''}{''}{''}{''}{'Acari'}{''}{''}{'Parasitiformes'}{''}{''}{''}{''}{''}{''}{''}{'Dermacentor'}{'Dermacentor_silvarum_(Dsil)'}{''} = 73;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{'Euchelicerata'}{''}{''}{''}{'Arachnida'}{''}{''}{''}{''}{''}{''}{'Arachnopulmonata'}{''}{''}{'Scorpiones'}{''}{''}{''}{''}{''}{''}{''}{'Centruroides'}{'Centruroides_sculpturatus_(Cscu)'}{''} = 74;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Chelicerata'}{'Euchelicerata'}{''}{''}{''}{'Arachnida'}{''}{''}{''}{''}{''}{''}{'Arachnopulmonata'}{''}{''}{'Araneae'}{''}{''}{''}{''}{''}{''}{''}{'Parasteatoda'}{'Parasteatoda_tepidariorum_(Ptep)'}{''} = 75;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Myriapoda'}{''}{''}{''}{'Diplopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Glomeris'}{'Glomeris_maerens_(Gmae)'}{''} = 76;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Myriapoda'}{''}{''}{''}{'Diplopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Helicorthomorpha'}{'Helicorthomorpha_holstii_(Hhol)'}{''} = 77;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Myriapoda'}{''}{''}{''}{'Chilopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Strigamia'}{'Strigamia_maritima_(Smar)'}{''} = 78;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Myriapoda'}{''}{''}{''}{'Chilopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Rhysida'}{'Rhysida_immarginata_(Rimm)'}{''} = 79;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Darwinula'}{'Darwinula_stevensoni_(Dste)'}{''} = 80;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Copepoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Eurytemora'}{'Eurytemora_affinis_(Eaff)'}{''} = 81;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Copepoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Tigriopus'}{'Tigriopus_californicus_(Tcal)'}{''} = 82;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Pollicipes'}{'Pollicipes_pollicipes_(Ppol)'}{''} = 83;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{'Malacostraca'}{''}{''}{''}{''}{''}{''}{''}{''}{'Decapoda'}{''}{''}{''}{''}{''}{''}{''}{'Penaeus'}{'Penaeus_chinensis_(Pchi)'}{''} = 84;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{'Malacostraca'}{''}{''}{''}{''}{''}{''}{''}{''}{'Decapoda'}{''}{''}{''}{''}{'Astacidea'}{''}{''}{'Homarus'}{'Homarus_americanusa_(Hame)'}{''} = 85;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{'Malacostraca'}{''}{''}{''}{''}{''}{''}{''}{''}{'Decapoda'}{''}{''}{''}{''}{'Astacidea'}{'Astacidea_a'}{''}{'Procambarus'}{'Procambarus_clarkii_(Pcla)'}{''} = 86;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{'Malacostraca'}{''}{''}{''}{''}{''}{''}{''}{''}{'Decapoda'}{''}{''}{''}{''}{'Astacidea'}{'Astacidea_a'}{''}{'Cherax'}{'Cherax_quadricarinatus_(Cqua)'}{''} = 87;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{'Malacostraca'}{''}{''}{''}{''}{''}{''}{''}{''}{'Multicrustacea_d'}{''}{''}{''}{''}{''}{''}{''}{'Armadillidium'}{'Armadillidium_nasatum_(Anas)'}{''} = 88;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Multicrustacea'}{'Multicrustacea_b'}{'Malacostraca'}{''}{''}{''}{''}{''}{''}{''}{''}{'Multicrustacea_d'}{''}{''}{''}{''}{''}{''}{''}{'Hyalella'}{'Hyalella_azteca_(Hazt)'}{''} = 89;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Arthropoda_d'}{'Branchiopoda'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Daphnia_m'}{'Daphnia_magna_(Dmag)'}{''} = 90;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Arthropoda_d'}{'Branchiopoda'}{'Branchiopoda_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Daphnia_g'}{'Daphnia_galeata_(Dgal)'}{''} = 91;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Arthropoda_d'}{'Branchiopoda'}{'Branchiopoda_a'}{'Branchiopoda_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Daphnia_x'}{'Daphnia_pulex_(Dpul)'}{''} = 92;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{''}{'Arthropoda_d'}{'Branchiopoda'}{'Branchiopoda_a'}{'Branchiopoda_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Daphnia_a'}{'Daphnia_pulicaria_(Dapu)'}{''} = 93;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{'Hexapoda'}{'Arthropoda_d'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Tomocerus'}{'Tomocerus_qinae_(Tqin)'}{''} = 94;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{'Hexapoda'}{'Arthropoda_d'}{'Hexapoda_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Campodea'}{'Campodea_augens_(Caug)'}{''} = 95;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{'Hexapoda'}{'Arthropoda_d'}{'Hexapoda_a'}{'Insecta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Ischnura'}{'Ischnura_elegans_(Iele)'}{''} = 96;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{'Hexapoda'}{'Arthropoda_d'}{'Hexapoda_a'}{'Insecta'}{'Insecta_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Zootermopsis'}{'Zootermopsis_nevadensis_(Znev)'}{''} = 97;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{'Hexapoda'}{'Arthropoda_d'}{'Hexapoda_a'}{'Insecta'}{'Insecta_a'}{'Insecta_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Tribolium'}{'Tribolium_castaneum_(Tcas)'}{''} = 98;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Protostomia'}{'Ecdysozoa'}{'Ecdysozoa_a'}{'Ecdysozoa_b'}{'Ecdysozoa_c'}{'Arthropoda'}{''}{'Arthropoda_a'}{'Arthropoda_b'}{'Arthropoda_c'}{'Hexapoda'}{'Arthropoda_d'}{'Hexapoda_a'}{'Insecta'}{'Insecta_a'}{'Insecta_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Drosophila'}{'Drosophila_melanogaster_(Dmel)'}{''} = 99;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Deuterostomia_a'}{'Hemichordata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Saccoglossus'}{'Saccoglossus_kowalevskii_(Skow)'}{''} = 100;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Deuterostomia_a'}{'Echinodermata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Anneissia'}{'Anneissia_japonica_(Anja)'}{''} = 101;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Deuterostomia_a'}{'Echinodermata'}{'Echinodermata_a'}{'Echinodermata_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Apostichopus'}{'Apostichopus_japonicus_(Ajap)'}{''} = 102;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Deuterostomia_a'}{'Echinodermata'}{'Echinodermata_a'}{'Echinodermata_b'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Strongylocentrotus'}{'Strongylocentrotus_purpuratus_(Spur)'}{''} = 103;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Deuterostomia_a'}{'Echinodermata'}{'Echinodermata_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Asterias'}{'Asterias_rubens_(Arub)'}{''} = 104;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Cephalochordata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Branchiostoma'}{'Branchiostoma_lanceolatum_(Blan)'}{''} = 105;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Tunicata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Ciona'}{'Ciona_intestinalis_(Cint)'}{''} = 106;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Myxini'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Eptatretus'}{'Eptatretus_burgeri_(Ebur)'}{''} = 107;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Chondrichthyes'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Callorhinchus'}{'Callorhinchus_milii_(Cmil)'}{''} = 108;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Chondrichthyes'}{'Chondrichthyes_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Leucoraja'}{'Leucoraja_erinacea_(Leri)'}{''} = 109;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Chondrichthyes'}{'Chondrichthyes_a'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Chiloscyllium'}{'Chiloscyllium_plagiosum_(Cpla)'}{''} = 110;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Danio'}{'Danio_rerio_(Drer)'}{''} = 111;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_c'}{''}{''}{''}{''}{''}{''}{'Sphaeramia'}{'Sphaeramia_orbicularis_(Sorb)'}{''} = 112;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_c'}{'Gobiiformes'}{''}{''}{''}{''}{''}{'Boleophthalmus'}{'Boleophthalmus_pectinirostris_(Bpec)'}{''} = 113;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_c'}{'Gobiiformes'}{'Gobiiformes_a'}{''}{''}{''}{''}{'Periophthalmus_a'}{'Periophthalmus_magnuspinnatus_(Pmag)'}{''} = 114;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_c'}{'Gobiiformes'}{'Gobiiformes_a'}{''}{''}{''}{''}{'Periophthalmus_o'}{'Periophthalmus_modestus_(Pmod)'}{''} = 115;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_b'}{'Actinopterygii_d'}{''}{''}{''}{''}{''}{'Oreochromis'}{'Oreochromis_niloticus_(Onil)'}{''} = 116;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_b'}{'Actinopterygii_d'}{''}{''}{''}{''}{''}{'Oryzias'}{'Oryzias_latipes_(Olat)'}{''} = 117;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Actinopterygii'}{''}{''}{''}{''}{''}{''}{''}{'Actinopterygii_a'}{'Actinopterygii_b'}{''}{''}{''}{''}{''}{''}{'Takifugu'}{'Takifugu_rubripes_(Trub)'}{''} = 118;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Actinistia'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Latimeria'}{'Latimeria_chalumnae_(Lcha)'}{''} = 119;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Dipnoi'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Protopterus'}{'Protopterus_annectens_(Pann)'}{''} = 120;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Gymnophiona'}{''}{''}{''}{''}{''}{''}{''}{'Rhinatrema'}{'Rhinatrema_bivittatum_(Rbiv)'}{''} = 121;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Gymnophiona'}{'Gymnophiona_a'}{''}{''}{''}{''}{''}{''}{'Geotrypetes'}{'Geotrypetes_seraphini_(Gser)'}{''} = 122;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Gymnophiona'}{'Gymnophiona_a'}{''}{''}{''}{''}{''}{''}{'Microcaecilia'}{'Microcaecilia_unicolor_(Muni)'}{''} = 123;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{''}{''}{''}{''}{''}{''}{''}{'Pleurodeles'}{'Pleurodeles_waltl_(Pwal)'}{''} = 124;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{''}{''}{''}{''}{''}{''}{'Bombina'}{'Bombina_bombina_(Bbom)'}{''} = 125;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{'Anura_a'}{''}{''}{''}{'Pipidae'}{''}{'Xenopus_l'}{'Xenopus_laevis_(Xlae)'}{''} = 126;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{'Anura_a'}{''}{''}{''}{'Pipidae'}{''}{'Xenopus_t'}{'Xenopus_tropicalis_(Xtro)'}{''} = 127;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{'Anura_a'}{'Anura_b'}{''}{''}{''}{''}{'Spea'}{'Spea_bombifrons_(Sbom)'}{''} = 128;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{'Anura_a'}{'Anura_b'}{''}{''}{'Anura_c'}{''}{'Bufo'}{'Bufo_bufo_(Bbuf)'}{''} = 129;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{'Anura_a'}{'Anura_b'}{''}{''}{'Anura_c'}{'Anura_d'}{'Nanorana'}{'Nanorana_parkeri_(Npar)'}{''} = 130;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Amphibia'}{''}{''}{''}{''}{'Amphibia_a'}{'Anura'}{'Anura_a'}{'Anura_b'}{''}{''}{'Anura_c'}{'Anura_d'}{'Rana'}{'Rana_temporaria_(Rtem)'}{''} = 131;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Squamata'}{''}{''}{''}{''}{''}{''}{''}{'Eublepharis'}{'Eublepharis_macularius_(Emac)'}{''} = 132;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Squamata'}{'Squamata_a'}{''}{''}{''}{''}{''}{''}{'Anolis'}{'Anolis_carolinensis_(Acar)'}{''} = 133;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Squamata'}{'Squamata_a'}{''}{''}{''}{''}{''}{''}{'Crotalus'}{'Crotalus_tigris_(Ctig)'}{''} = 134;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{''}{''}{''}{''}{''}{''}{'Pelodiscus'}{'Pelodiscus_sinensis_(Psin)'}{''} = 135;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_b'}{''}{''}{'Testudines_d'}{''}{'Dermochelys'}{'Dermochelys_coriacea_(Deco)'}{''} = 136;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_b'}{''}{''}{'Testudines_d'}{'Cheloniidae'}{'Caretta'}{'Caretta_caretta_(Care)'}{''} = 137;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_b'}{''}{''}{'Testudines_d'}{'Cheloniidae'}{'Chelonia'}{'Chelonia_mydas_(Cmyd)'}{''} = 138;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_b'}{''}{''}{''}{''}{'Chelydra'}{'Chelydra_serpentina_(Cser)'}{''} = 139;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_c'}{''}{''}{''}{''}{'Platysternon'}{'Platysternon_megacephalum_(Pmeg)'}{''} = 140;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_c'}{''}{''}{'Emydidae'}{''}{'Chrysemys'}{'Chrysemys_picta_bellii_(Cpic)'}{''} = 141;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Testudines'}{'Testudines_a'}{'Testudines_c'}{''}{''}{'Emydidae'}{''}{'Terrapene'}{'Terrapene_carolina_triunguis_(Tcar)'}{''} = 142;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Vertebrata_i'}{'Crocodylia'}{''}{''}{''}{'Alligatoridae'}{''}{'Alligator_m'}{'Alligator_mississippiensis_(Amis)'}{''} = 143;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Vertebrata_i'}{'Crocodylia'}{''}{''}{''}{'Alligatoridae'}{''}{'Alligator_s'}{'Alligator_sinensis_(Asin)'}{''} = 144;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Vertebrata_i'}{'Crocodylia'}{'Crocodylia_a'}{''}{''}{''}{''}{'Crocodylus'}{'Crocodylus_porosus_(Cpor)'}{''} = 145;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Vertebrata_i'}{'Crocodylia'}{'Crocodylia_a'}{''}{''}{''}{''}{'Gavialis'}{'Gavialis_gangeticus_(Ggan)'}{''} = 146;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Vertebrata_i'}{'Aves'}{''}{''}{''}{''}{''}{'Gallus'}{'Gallus_gallus_(Ggal)'}{''} = 147;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Vertebrata_g'}{''}{''}{''}{'Vertebrata_h'}{'Vertebrata_i'}{'Aves'}{''}{''}{''}{''}{''}{'Taeniopygia'}{'Taeniopygia_guttata_(Tgut)'}{''} = 148;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Mammalia'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Ornithorhynchus'}{'Ornithorhynchus_anatinus_(Oana)'}{''} = 149;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Mammalia'}{''}{''}{''}{'Mammalia_a'}{''}{''}{''}{''}{''}{''}{''}{'Monodelphis'}{'Monodelphis_domestica_(Mdom)'}{''} = 150;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Mammalia'}{''}{''}{''}{'Mammalia_a'}{'Mammalia_b'}{''}{''}{''}{''}{''}{''}{'Homo'}{'Homo_sapiens_(Hsap)'}{''} = 151;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Mammalia'}{''}{''}{''}{'Mammalia_a'}{'Mammalia_b'}{'Mammalia_c'}{''}{''}{''}{''}{''}{'Canis'}{'Canis_lupus_familiaris_(Clup)'}{''} = 152;
$spp {'Eukaryota'}{'Eukaryota_a'}{'Eukaryota_b'}{'Metazoa'}{'Metazoa_a'}{'Planulozoa'}{'Bilateria'}{'Nephrozoa'}{'Deuterostomia'}{'Chordata'}{''}{''}{''}{''}{''}{'Chordata_a'}{'Vertebrata'}{''}{''}{''}{'Vertebrata_a'}{'Vertebrata_b'}{'Vertebrata_c'}{'Vertebrata_d'}{'Vertebrata_e'}{'Vertebrata_f'}{'Mammalia'}{''}{''}{''}{'Mammalia_a'}{'Mammalia_b'}{'Mammalia_c'}{''}{''}{''}{''}{''}{'Tursiops'}{'Tursiops_truncatus_(Ttru)'}{''} = 153;

print_hash_colors (%spp);
return %spp;
}

#PRINT_HASH_COLORS: subroutine to traverse the hash of hashes and print the contents, after http://www.perlmonks.org/?node_id=116162
sub print_hash_colors {
    my (%spp) = @_; 
#Define variables (taxonomic ranks)for the hash of hashes
my $domain;
my $node_a;
my $node_b;
my $node_c;
my $phylum_a;
my $phylum_b;
my $phylum_c;
my $phylum_d;
my $phylum_e;
my $phylum_f;
my $phylum_g;
my $phylum_h;
my $phylum_i;
my $phylum_j;
my $phylum_k;
my $subphylum_a;
my $subphylum_b;
my $subphylum_c;
my $subphylum_d;
my $superclass;
my $class_a;
my $class_b;
my $class_c;
my $class_d;
my $class_e;
my $class_f;
my $class_g;
my $subclass_a;
my $subclass_b;
my $subclass_c;
my $order_a;
my $order_b;
my $order_c;
my $order_d;
my $suborder;
my $infraorder;
my $family_a;
my $family_b;
my $genus;
my $species;

#Print all the keys of the hash of hashes with a nested FOR loop
###   domain # node_a # node_b # node_c # phylum_a # phylum_b # phylum_c # phylum_d # phylum_e # phylum_f # phylum_g # phylum_h # phylum_i # phylum_j # phylum_k # subphylum_a # subphylum_b # subphylum_c # subphylum_d # superclass # class_a # class_b # class_c # class_d # class_e # class_f # class_g # subclass_a # subclass_b # subclass_c # order_a # order_b # order_c # order_d # suborder # infraorder # family_a # family_b # genus # species

print "\nTree:";
for $domain (keys %spp) {
    print "\n$domain\n";
    for $node_a (keys %{ $spp{$domain} }) {
        unless ($node_a eq '') {print color ("blue"),"$node_a\n", color ("reset");}
        for $node_b (keys %{ $spp{$domain}{$node_a} }) {
            unless ($node_b eq '') { print color ("yellow"),"  $node_b\n", color ("reset");}
            for $node_c (keys %{ $spp{$domain}{$node_a}{$node_b} }) {
                unless ($node_c eq '') { print color ("green"),"    $node_c\n", color ("reset");}
                for $phylum_a (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c} }) {
                    unless ($phylum_a eq '') { print color ("magenta"),"      $phylum_a\n", color ("reset");}
                    for $phylum_b (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a} }) {
                        unless ($phylum_b eq '') { print color ("bright_cyan"),"        $phylum_b\n", color ("reset");}
                        for $phylum_c (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b} }) {
                            unless ($phylum_c  eq '') { print color ("red"),"           $phylum_c\n", color ("reset");}
                            for $phylum_d (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c} }) {
                                unless ($phylum_d eq '') { print color ("blue"),"             phylum_d\n", color ("reset");}
                                for $phylum_e (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d} }) {
                                    unless ($phylum_e eq '') { print color ("bright_yellow"),"               phylum_e\n", color ("reset");}
                                    for $phylum_f (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e} }) {
                                        unless ($phylum_f eq '') { print color ("bright_green"),"                 phylum_f\n", color ("reset");}
                                        for $phylum_g (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f} }) {
                                            unless ($phylum_g eq '') { print color ("bright_magenta"),"                   phylum_g\n", color ("reset");}
                                            for $phylum_h (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g} }) {
                                                unless ($phylum_h eq '') { print color ("bright_cyan"),"                     phylum_h\n", color ("reset");}
                                                for $phylum_i (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h} }) {
                                                    unless ($phylum_i eq '') { print color ("bright_red"),"                       phylum_i\n", color ("reset");}
                                                    for $phylum_j (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i} }) {
                                                        unless ($phylum_j eq '') { print color ("bright_blue"),"                         phylum_j\n", color ("reset");}
                                                        for $phylum_k (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j} }) {
                                                            unless ($phylum_k eq '') { print color ("blue"),"                           phylum_k\n", color ("reset");}
                                                            for $subphylum_a (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k} }) {
                                                                unless ($subphylum_a eq '') { print color ("yellow"),"                             subphylum_a\n", color ("reset");}
                                                                for $subphylum_b (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a} }) {
                                                                    unless ($subphylum_b eq '') { print color ("green"),"                               subphylum_b\n", color ("reset");}
                                                                    for $subphylum_c (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b} }) {
                                                                        unless ($subphylum_c eq '') { print color ("magenta"),"                               subphylum_c\n", color ("reset");}
                                                                        for $subphylum_d (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c} }) {
                                                                            unless ($subphylum_d eq '') { print color ("bright_cyan"),"                                 subphylum_d\n", color ("reset");}
                                                                            for $superclass (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d} }) {
                                                                                unless ($superclass eq '') { print color ("red"),"                                   superclass\n", color ("reset");}
                                                                                for $class_a (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass} }) {
                                                                                    unless ($class_a eq '') { print color ("blue"),"                                     class_a\n", color ("reset");}
                                                                                    for $class_b (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a} }) {
                                                                                        unless ($class_b eq '') { print color ("bright_yellow"),"                                       class_b\n", color ("reset");}
                                                                                        for $class_c (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b} }) {
                                                                                            unless ($class_c eq '') { print color ("bright_green"),"                                         class_c\n", color ("reset");}
                                                                                            for $class_d (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c} }) {
                                                                                                unless ($class_d eq '') { print color ("bright_magenta"),"                                           class_d\n", color ("reset");}
                                                                                                for $class_e (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d} }) {
                                                                                                    unless ($class_e eq '') { print color ("bright_cyan"),"                                             class_e\n", color ("reset");}
                                                                                                    for $class_f (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e} }) {
                                                                                                        unless ($class_f eq '') { print color ("bright_red"),"                                               class_f\n", color ("reset");}
                                                                                                        for $class_g (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f} }) {
                                                                                                            unless ($class_g eq '') { print color ("bright_blue"),"                                                 class_g\n", color ("reset");}
                                                                                                            for $subclass_a (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g} }) {
                                                                                                                unless ($subclass_a eq '') { print color ("blue"),"                                                   subclass_a\n", color ("reset");}
                                                                                                                for $subclass_b (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a} }) {
                                                                                                                    unless ($subclass_b eq '') { print color ("yellow"),"                                                   subclass_b\n", color ("reset");}
                                                                                                                    for $subclass_c (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b} }) {
                                                                                                                        unless ($subclass_c eq '') { print color ("green"),"                                                     subclass_c\n", color ("reset");}
                                                                                                                        for $order_a (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c} }) {
                                                                                                                            unless ($order_a eq '') { print color ("magenta"),"                                                       order_a\n", color ("reset");}
                                                                                                                            for $order_b (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a} }) {
                                                                                                                                unless ($order_b eq '') { print color ("bright_cyan"),"                                                         order_b\n", color ("reset");}
                                                                                                                                for $order_c (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b} }) {
                                                                                                                                    unless ($order_c eq '') { print color ("red"),"                                                           order_c\n", color ("reset");}
                                                                                                                                    for $order_d (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c} }) {
                                                                                                                                        unless ($order_d eq '') { print color ("blue"),"                                                             order_d\n", color ("reset");}
                                                                                                                                        for $suborder (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d} }) {
                                                                                                                                            unless ($suborder eq '') { print color ("bright_yellow"),"                                                               suborder\n", color ("reset");}
                                                                                                                                            for $infraorder (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d}{$suborder} }) {
                                                                                                                                                unless ($infraorder eq '') { print color ("bright_green"),"                                                                 infraorder\n", color ("reset");}
                                                                                                                                                for $family_a (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d}{$suborder}{$infraorder} }) {
                                                                                                                                                    unless ($family_a eq '') { print color ("bright_magenta"),"                                                                   family_a\n", color ("reset");}
                                                                                                                                                    for $family_b (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d}{$suborder}{$infraorder}{$family_a} }) {
                                                                                                                                                        unless ($family_b eq '') { print color ("bright_cyan"),"                                                                     family_b\n", color ("reset");}
                                                                                                                                                        for $genus (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d}{$suborder}{$infraorder}{$family_a}{$family_b} }) {
                                                                                                                                                            unless ($genus eq '') { print color ("bright_red"),"                                                                       genus\n", color ("reset");}
                                                                                                                                                            for $species (keys %{ $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d}{$suborder}{$infraorder}{$family_a}{$family_b}{$genus} }) {
                                                                                                                                                                print color ("bright_white"),"                                                                         $species => $spp{$domain}{$node_a}{$node_b}{$node_c}{$phylum_a}{$phylum_b}{$phylum_c}{$phylum_d}{$phylum_e}{$phylum_f}{$phylum_g}{$phylum_h}{$phylum_i}{$phylum_j}{$phylum_k}{$subphylum_a}{$subphylum_b}{$subphylum_c}{$subphylum_d}{$superclass}{$class_a}{$class_b}{$class_c}{$class_d}{$class_e}{$class_f}{$class_g}{$subclass_a}{$subclass_b}{$subclass_c}{$order_a}{$order_b}{$order_c}{$order_d}{$suborder}{$infraorder}{$family_a}{$family_b}{$genus}{$species}{''}\n", color ("reset");
                                                                                                                                                            }
                                                                                                                                                        }
                                                                                                                                                    }
                                                                                                                                                }
                                                                                                                                            }
                                                                                                                                        }
                                                                                                                                    }
                                                                                                                                }
                                                                                                                            }
                                                                                                                        }
                                                                                                                    }
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}                    
return;
}