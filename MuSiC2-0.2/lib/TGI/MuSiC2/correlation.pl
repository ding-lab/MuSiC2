#!/opt/local/bin/perl

#  FOR EACH GENE, DETERMINE THE CORRELATION BETWEEN ITS CANCER-SPECIFIC
#  -LOG10(P) AND THAT CANCER'S RESCALED BACKGROUND MUTATION RATE AND ALSO THE
#  MAXIMUM DISTANCE OF POINTS TO THE REGRESSION LINE
#
#  REFERENCE completed/2007_correction_factor_fpc/scripts/plot_data.pl FOR SOME
#  PROGRAMMING GOODIES RELATIVE TO THIS PROBLEM
#
#  IN READING THE DATA, THE LINE "if (defined $pval_neg_log && $pval_neg_log) {"
#  **IS** READING OCCURENCES OF 0.0 IN THE INPUT FILES BECAUSE, ALTHOUGH PERL
#  PROCESSES NUMBERS AND TEXT-THAT-RESEMBLES-A-NUMBER DIFFERENTLY (SEE BELOW
#  EXAMPLE), IT *READS* A FILE OF NUMBERS INITIALLY AS TEXT
#
#   bash-3.2$ perl             bash-3.2$ perl
#   $s = 0.0;                  $s = "0.0";   <---quotation marks on this one
#   if ($s) {                  if ($s) {
#     print "yes\n";              print "yes\n";
#   } else {                   } else {
#     print "no\n";               print "no\n";
#   }                          }
#   no                         yes
#
#   USAGE> ./correlation.pl > plot_file.dat

#  INFER CANCER TYPE FROM P-VALUE FILENAME AND VERIFY CONSISTENCY WITH RATE
#  NOMENCLATURE --- THIS DOES NOT INFER THE CANCER TYPE OF THE GENES-OF-INTEREST
#  FILE --- THAT IS SEPARATE
# http://stackoverflow.com/questions/20647010/pass-regex-into-perl-subroutine
# my $regexp = qr/^Perl$/;

############
#  SET UP  #
############

#__STANDARD PERL PACKAGES
   use strict;
   use warnings;

#__STATISTICAL LIBRARIES
   use Statistics::OLS;

###################
#  CONFIGURATION  #
###################

#__FILE OF GENE NAMES FOR WHICH TO RUN ANALYSIS
#
#  PROGRAMMING NOTE: SOME FILES FROM BAILEY (FROM MICROSOFT APPS) HAVE
#  HAVE BOTH \N AND \CR SO **ALWAYS** CHECK THE READING LOOP BELOW TO DETERMINE
#  WHETHER "chop" IS NEEED
#
#  CANCER TYPE IS INFERRED FROM THE FILENAME: SHOULD HAVE FORM OF "CANCER_*"

#  my $file_gene_list = "laml.tophits";
#  my $file_gene_list = "321_Round/STAD_SMGs_from_321_genes.txt";
#  my $file_gene_list = "LUAD_SMGs_from_366_genes_with_max_can.txt";
   my $file_gene_list = "plots/cc_vs_cancer_dis_all_genes_annotated_from_366/STAD_SMGs_from_366_genes.txt";

#__ANALYZE "KANDOTH" SUBSET OF THESE GENES (=1) OR THE "NOT KANDOTH" SUBSET (=0)
   my $in_kandoth = 0;

#__PANCAN DATA FILES TO PROCESS
   my @files = qw(
      ../revised_manuscript_adjust_cutoff/acc_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/blca_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/brca_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/cesc_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/coadread_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/gbm_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/hnsc_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/kich_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/kirc_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/kirp_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/laml_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/lgg_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/lihc_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/luad_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/lusc_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/OV_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/prad_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/skcm_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/stad_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/thca_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/ucec_pval.genesize.txt
      ../revised_manuscript_adjust_cutoff/ucs_pval.genesize.txt
   );

#__CANCER BACKGROUND MUTATION RATES (TRANSFORM TO 25-SCALE BY MULT 8.5e6)
   my $cancer_background_rate = {
      'acc' =>      5.04443724040595e-07 * 8500000,
      'blca' =>     2.29927147645268e-06 * 8500000,
      'brca' =>     4.43871323686807e-07 * 8500000,
      'cesc' =>     9.58381886685746e-07 * 8500000,
      'coadread' => 8.1086390892704e-07 * 8500000,
      'gbm' =>      5.59470325737927e-07 * 8500000,
      'hnsc' =>     1.29729704730899e-06 * 8500000,
      'kich' =>     6.43549698064298e-07 * 8500000,
      'kirc' =>     5.47376310309256e-07 * 8500000,
      'kirp' =>     9.19530214627311e-07 * 8500000,
      'laml' =>     1.16867075954987e-07 * 8500000,
      'lgg' =>      3.86862498457639e-07 * 8500000,
      'lihc' =>     1.95329184153507e-06 * 8500000,
      'luad' =>     2.3821583999372e-06 * 8500000,
      'lusc' =>     2.68752209654077e-06 * 8500000,
      'ov' =>       5.94633670211403e-07 * 8500000,
      'prad' =>     3.62495550231542e-07 * 8500000,
      'skcm' =>     2.90983510858331e-06 * 8500000,
      'stad' =>     1.52211643106822e-06 * 8500000,
      'thca' =>     1.46669747383227e-07 * 8500000,
      'ucec' =>     1.54617002760412e-06 * 8500000,
      'ucs' =>      6.72611552630032e-07 * 8500000,
   };

####################
#  PRE-PROCESSING  #
####################

#__DETERMINE OUTPUT FILE NAME
   my $file_out;
   if ($file_gene_list =~ /^(\S+)\.txt$/) {
      $file_out = $1;
   } else {
      die "unexpected naming for input file '$file_gene_list'";
   }
   if ($in_kandoth) {
      $file_out .= "_kandoth_YES.dat";
   } else {
      $file_out .= "_kandoth_NO.dat";
   }
   print "OUTPUT FILE IS  $file_out\n";

#__OPEN OUTPUT FILE
   unlink $file_out if -e $file_out;
   open (OUT, ">$file_out") || die "cant open $file_out";

#__READ P-VALUE DATA FOR EACH GENE ACROSS ALL AVAILABLE CANCER TYPES
   my ($gene_pval_data, $gene_size) = ({}, {});
   foreach my $file (@files) {

   #__INFER CANCER TYPE FROM FILE NAME USING REGEXP (SEE PROGRAMMING NOTES)
      my $regexp = qr/(\w+)\_pval\.genesize\.txt$/;
      my $cancer = _cancer_type_ ($file, $regexp);

   #__READ DATA
      open (F, $file) || die "cant open $file";
      while (<F>) {

      #__PARSE PROCESS AND STORE A LINE OF DATA
         next if /^#/;
         chomp;
         my ($gene, $size, $pval_neg_log) = split;
         if (defined $pval_neg_log && $pval_neg_log) {
            $gene_pval_data->{$gene}->{$cancer} = $pval_neg_log;
            $gene_size->{$gene} = $size; # OVER-WRITTEN WITH SAME DATA
         }
      }
      close (F);
   }

#__READ KANDOTH LIST OF BONA-FIDE CANCER GENES
   my $kandoth = {};
   while (<DATA>) {
      chomp;
      $kandoth->{$_} = 1;
   }

#__READ THE LIST OF GENES OF INTEREST
   my $gene_list = {};
   open (F, $file_gene_list) || die "cant open $file_gene_list";
   while (<F>) {
      chomp;
#     chop; # 943genesPancan2_from_Bailey.txt SEEMS TO HAVE BOTH \N AND \CR
      $gene_list->{$_} = 1;
#     warn "reading --->$_<---\n";
   }
   close (F);

#__INFER CANCER TYPE FROM GENES-OF-INTEREST FILE USING REGEXP (SEE PROG NOTES)
   my $regexp = qr/^plots\/cc_vs_cancer_dis_all_genes_annotated_from_366\/([A-Za-z]+)\_/;
   my $cancer_type_main = _cancer_type_ ($file_gene_list, $regexp);

#####
#  close (OUT);
#  die "cancer is '$cancer_type_main' --- END";
#####

#__HEADER
   print OUT "# CORRELATION COEFFICIENT VS CANCER-SPECIFIC DISTANCE TO REGRESSION LINE\n#\n";
   print OUT "# script: $0\n#\n";
   print OUT "# gene list input file: $file_gene_list\n#\n";
   print OUT "# cancer type: $cancer_type_main\n#\n";
   if ($in_kandoth) {
      print OUT "# status of genes being screened: in the Kandoth subset\n#\n";
   } else {
      print OUT "# status of genes being screened: NOT in the Kandoth subset\n#\n";
   }

#####################
#  MAIN PROCESSING  #
#####################
#
#  X-AXIS: CANCER BACKGROUND MUTATION RATE * 8.5e+6
#  Y-AXIS: -LOG10(P-VALUE)

#__GENE-BY-GENE REGRESSION ANALYSIS
   my ($x_axis_min, $x_axis_max) = (9999, 0);
   my ($y_axis_min, $y_axis_max) = (9999, 0);
   my $num_filtered = 0;
   foreach my $gene (sort keys %{$gene_pval_data}) {

   #__WE ARE ONLY LOOKING AT GENES IN THE LIST OF INTEREST
      next unless defined $gene_list->{$gene};

   #__GENERATE THE INPUT XVALS AND YVALS VECTORS
      my ($num_non_0_yvals, $num_max_significant_p) = (0, 0);
      my ($xvals, $yvals) = ([], []);
      foreach my $cancer (sort keys %{$gene_pval_data->{$gene}}) {
         $num_non_0_yvals++ if $gene_pval_data->{$gene}->{$cancer} > 0;
         $num_max_significant_p++ if $gene_pval_data->{$gene}->{$cancer} > 15;
         push @{$xvals}, $cancer_background_rate->{$cancer};
         push @{$yvals}, $gene_pval_data->{$gene}->{$cancer};
      }
      my $num_pts = scalar @{$xvals};

   #__MAKE SURE GENE HAS AT LEAST 3 DATA POINTS
      if ($num_pts < 3) {
         warn "$gene HAS LESS THAN 3 POINTS --- SKIPPING\n";
         next;
      }

   #__EXECUTE REGRESSION
      my $stats_obj = Statistics::OLS->new();
      $stats_obj->setData ($xvals, $yvals) || die ($stats_obj->error());
      $stats_obj->regress () || die ($stats_obj->error());

   #__PEARSON'S CORRELATION COEFFICIENT INDICATES GOODNESS OF CORRELATION
      my $r_squared = $stats_obj->rsq();
      $r_squared = 0 if $r_squared < 0;
      my $pearson = sqrt $r_squared;

   #__CALCULATE REGRESSION EQUATION PARAMETERS
      my ($intercept, $slope) = $stats_obj->coefficients();
      if ($slope == 0) {
         warn "$gene HAS 0 SLOPE --- CHECK YOUR DATA --- SKIPPING\n";
         next;
      }

   #__FIND MAX DISTANCE OF A POINT TO THE REGRESSION LINE (SEE OFFLINE NOTES)
      my ($y_value, $y_cancer) = (0, "");
      foreach my $cancer (sort keys %{$gene_pval_data->{$gene}}) {
         my ($x_1, $y_1) = ($cancer_background_rate->{$cancer}, $gene_pval_data->{$gene}->{$cancer});
         my $x_2 = ($y_1 - $intercept + $x_1/$slope) / ($slope + 1/$slope);
         my $y_2 = $slope*$x_2 + $intercept;
         my $distance = sqrt(($x_2-$x_1)*($x_2-$x_1) + ($y_2-$y_1)*($y_2-$y_1));

      ##_Y-AXIS IS MAXIMUM DISTANCE TO REGRESSION LINE (DEPRECATED)
      ## if ($distance > $y_value) {
      ##    $y_value = $distance;
      ##    $y_cancer = $cancer;
      ## }

      #__Y-AXIS IS DISTANCE TO REGRESSION LINE FOR CANCER-SPECIFIC DISTANCE
         if ($cancer eq $cancer_type_main) {
            $y_value = $distance;
            $y_cancer = $cancer; # TECHNICALLY SUPERFLUOUS
         }
      }

   #__CONDITIONAL OUTPUT
      if ($in_kandoth) {

      #__FOR KANDOTH GENES WE MAY OR MAY NOT NEED GENE NAMES ANNOTATED
         if (defined $kandoth->{$gene}) {
      ##    print OUT "$pearson    $y_value\n";
            print OUT "$pearson    $y_value   \"$gene\"\n";
      ####  print OUT "$pearson    $y_value   \"$gene   $y_cancer\"\n";
         }
      } else {

      #__FOR NON-KANDOTH GENES WE NEED GENE NAMES ANNOTATED
         unless (defined $kandoth->{$gene}) {
            print OUT "$pearson    $y_value   \"$gene\"\n";
      ####  print OUT "$pearson    $y_value   \"$gene   $y_cancer\"\n";
         }
      }

   #__TRACK AXIS BOUNDARIES TO REPORT HOW BIG X-AXIS AND Y-AXIS SHOULD BE SET TO
      $x_axis_max = $pearson if $pearson > $x_axis_max;
      $x_axis_min = $pearson if $pearson < $x_axis_min;
      $y_axis_max = $y_value if $y_value > $y_axis_max;
      $y_axis_min = $y_value if $y_value < $y_axis_min;
   }

#__OUTPUT AXIS BOUNDARIES FOR PLOTTING: COMMENTS IGNORED BY XMGRACE
   print OUT "#  $x_axis_min <= X <= $x_axis_max\n";
   print OUT "#  $y_axis_min <= Y <= $y_axis_max\n";

##################
#  POST PROCESS  #
##################

#__CLOSE OUTPUT FILE
   close (OUT);

################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################

#  INFER CANCER TYPE EMBEDDED IN A FILE NAME: USES COMPILED REGEXPS HAVING
#  A SINGLE CALLBACK --- SEE PROGRAMMING NOTES ABOVE

sub _cancer_type_ {
   my ($file, $regexp) = @_;
   my $cancer;
   if ($file =~ /$regexp/) {
      $cancer = $1;

   #__CONVERT MIXED TEXT TO ALL LOWER CASE FOR CONSISTENCY
      $cancer = lc $cancer;

   #__CHECK CONSISTENCY WITH BACKGROUND MUTATION RATE NOMENCLATURE
      unless (defined $cancer_background_rate->{$cancer}) {
        die "cancer '$cancer' from filename '$file' not a valid cancer type";
      }
   } else {
      die "cannot infer cancer type from file name '$file'";
   }
   return $cancer;
}

################################################################################
#                                                                              #
#                                   D A T A                                    #
#                                                                              #
################################################################################

#  THE GENES FROM KANDOTH ET AL THAT ARE BONA FIDE CANCER GENES --- THESE
#  ARE USEFUL TO GET A FEEL FOR WHAT PEARSON COEFFICIENT AND SLOPE ARE
#  CHARACTERISTIC OF REAL CANCER GENES --- THIS LIST IS **IDENTICAL** TO
#  THE CONTENTS OF FILE
#
#           Kandoth_127_genes_from_Mike_M.txt
#
#  BECAUSE WE READ-IN THE CONTENTS FROM THAT FILE

__DATA__
ACVR1B
ACVR2A
AJUBA
AKT1
APC
AR
ARHGAP35
ARID1A
ARID5B
ASXL1
ATM
ATR
ATRX
AXIN2
B4GALT3
BAP1
BRAF
BRCA1
BRCA2
CBFB
CCND1
CDH1
CDK12
CDKN1A
CDKN1B
CDKN2A
CDKN2C
CEBPA
CHEK2
CRIPAK
CTCF
CTNNB1
DNMT3A
EGFR
EGR3
EIF4A2
ELF3
EP300
EPHA3
EPHB6
EPPK1
ERBB4
ERCC2
EZH2
FBXW7
FGFR2
FGFR3
FLT3
FOXA1
FOXA2
GATA3
H3F3C
HGF
HIST1H1C
HIST1H2BD
IDH1
IDH2
KDM5C
KDM6A
KEAP1
KIT
KMT2B
KMT2C
KMT2D
KRAS
LIFR
LRRK2
MALAT1
MAP2K4
MAP3K1
MAPK8IP1
MECOM
MIR142
MTOR
NAV3
NCOR1
NF1
NFE2L2
NFE2L3
NOTCH1
NPM1
NRAS
NSD1
PBRM1
PCBP1
PDGFRA
PHF6
PIK3CA
PIK3CG
PIK3R1
POLQ
PPP2R1A
PRX
PTEN
PTPN11
RAD21
RB1
RPL22
RPL5
RUNX1
SETBP1
SETD2
SF3B1
SIN3A
SMAD2
SMAD4
SMC1A
SMC3
SOX17
SOX9
SPOP
STAG2
STK11
TAF1
TBL1XR1
TBX3
TET2
TGFBR2
TLR4
TP53
TSHZ2
TSHZ3
U2AF1
USP9X
VEZF1
VHL
WT1
