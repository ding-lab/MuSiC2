#!/opt/local/bin/perl

#  FIND CONDITIONS FOR WHICH SIGNIFICANCE STATUS IS NO LONGER RELATED TO GENE
#  SIZE --- THIS IS DONE BY ITERATIVELY RAISING THE CUT-OFF FOR LONGER GENES
#  AS COMPARED TO SHORTER ONES
#
#  =================
#  PROGRAMMING NOTES
#  =================
#
#  (1) GOODNESS OF FIT 2X2 TABLE TEST FOR EACH CANCER IN THE FORM OF
#      OF THE NUMBERS OF GENES IN EACH OF 2 CATEGORIES: (LONG OR NOT LONG) AND
#      (SIGNIFICANT P-VALUE OR NOT SIGNIFICANT P-VALUE)
#
#                          | NOT LONG     LONG |
#       -------------------+-------------------+---------------
#                          |                   |
#       SIGNIF P-VALUE     |   a          b    |   a + b
#                          |                   |
#       NOT SIGNIF P-VALUE |   c          d    |   c + d
#                          |                   |
#       -------------------+-------------------+---------------
#                          | a + c      b + d  | a + b + c + d
#
#   NOTES:
#
#   * note the protocal for using Statistics::Distributions as written in
#     research/in_progress/integrated_analysis/scripts/Statistics/Normality.pm
#
#   * from Statistics::Distributions internal documentation
#
#     $chisprob=Statistics::Distributions::chisqrprob (3,6.25);
#     print "upper probability of the chi-square distribution (3 degrees "
#            ."of freedom, chi-squared = 6.25): Q = 1-G = $chisprob\n";

############
#  SET UP  #
############

#__STANDARD PERL PACKAGES
   use strict;
   use warnings;

#__LIBRARIES
   use Statistics::Distributions;

###################
#  CONFIGURATION  #
###################

#__PANCAN DATA FILES TO PROCESS
#  my @files = qw/luad_pval.genesize.txt/;
   my @files = qw/OV_pval.genesize.txt   acc_pval.genesize.txt
                  blca_pval.genesize.txt brca_pval.genesize.txt
                  cesc_pval.genesize.txt coadread_pval.genesize.txt
                  gbm_pval.genesize.txt  hnsc_pval.genesize.txt
                  kich_pval.genesize.txt kirc_pval.genesize.txt
                  kirp_pval.genesize.txt laml_pval.genesize.txt
                  lgg_pval.genesize.txt  lihc_pval.genesize.txt
                  luad_pval.genesize.txt lusc_pval.genesize.txt
                  prad_pval.genesize.txt skcm_pval.genesize.txt
                  stad_pval.genesize.txt thca_pval.genesize.txt
                  ucec_pval.genesize.txt ucs_pval.genesize.txt/;

#__SET THE THRESHOLD BETWEEN "TYPICAL" GENES AND "LARGE" GENES
   my $x_thresh_gene_size = 5000;
#  my $x_thresh_gene_size = 4000;

#__DOWN-INCREMENT OF THRESHOLD DELIMITING TYPICAL AND LARGE GENES
   my $x_incr = 500;

#__HIGHEST P-VALUE FOR TYPICAL GENES IS ALSO LOWER-BOUND FOR LARGE GENES
#  (THE ACTUAL P-VALUE IS EXP (- $y_thresh_fixed)
#  my $y_thresh_fixed = 5; # 0.00001
#  my $y_thresh_fixed = 6; # 0.000001
#  my $y_thresh_fixed = 7; # 0.0000001
   my $y_thresh_fixed = 8; # 0.00000001

#__SET ABSOLUTE LOWEST ALLOWABLE P-VALUE (HIGHEST Y-THRESHOLD) TO USE
#  my $y_max = 13.816;
#  my $y_max = 16.118;
   my $y_max = 18;

#__INCREMENT OF LN(P-VALUE)
   my $y_incr = 0.1;

#__P-VALUE THRESHOLD ABOVE WHICH WE SAY THAT MUTATIONAL SIGNIFICANCE STATUS
#  IS NO LONGER RELATED TO GENE SIZE
#  my $pvalue_threshold = 0.01;
   my $pvalue_threshold = 0.005;
#  my $pvalue_threshold = 0.001;
#  my $pvalue_threshold = 0.0065;
#  my $pvalue_threshold = 0.05;

####################
#  PRE-PROCESSING  #
####################

#__HEADER
   print "# CANCER-SPECIFIC ADJUSTMENT OF P-VAL THRESHOLD FOR LARGE GENES\n#\n";
   print "# script: $0\n#\n";
   print "# ************************************************************\n#\n";
   print "# PARAMETERS\n#\n";
   print "# default delineation between 'typical' and 'large' genes: $x_thresh_gene_size\n#\n";
   print "#    (note that this boundary may be modified for some cancers)\n#\n";
   print "# highest allowable Y-threshold: $y_max\n#\n";
   print "# P-value target above which mutational significance status\n";
   print "#      and gene size are no longer statistically related: $pvalue_threshold\n#\n";
   print "# ************************************************************\n#\n";

#__BONRFERRONI CORRECTION BASED ON NUMBER OF CANCERS (DEPRECATED)
   my $num_cancers = scalar @files;
#  $pvalue_threshold /= $num_cancers;

#__GET LIST OF POSITIVE GENES (FOR CALLER ASSESSMENT)
   my $pos_genes = {};
   while (<DATA>) {
      chomp;
      $pos_genes->{$_} = 1;
   }

#####################
#  MAIN PROCESSING  #
#####################

#__PROCESS EACH FILE
   my $all_filtered_long_genes = {};
   my $total_gene_evals = 0;
   foreach my $file (@files) {
      my $eval_done = 0;

   #__READ DATA IN THIS FILE
      warn "# processing data in file: $file\n";
#     print "# -------------------------------------------------------------\n";
      print "# processing data in file: $file\n";
      print "# -------------------------------------------------------------\n";
      my ($pts, $gene_data, $max_y_large_gene) = _read_file_ ($file);

   #__Y-THRESHOLD IS CURRENTLY FIXED FOR LEFT-PLANE GENES
   #  my $y_thresh_l = $y_thresh_fixed;
      my $iters = 0;
      for (my $y_thresh_l = $y_thresh_fixed; $y_thresh_l >= 6; $y_thresh_l -= 0.5) {
         warn "#   trying cutoff of Y = $y_thresh_l\n";

   #__X-THRESHOLD THAT DELINEATES TYPICAL VS LARGE GENES
   #  my ($iters, $convergence) = (0, 0);
      for (my $x_thresh = $x_thresh_gene_size; $x_thresh > $x_incr;
                                               $x_thresh -= $x_incr) {

      #__MOVE Y-THRESHOLD UPWARD FOR LARGE GENES UNTIL SIGNIFICANCE VANISHES
         for (my $y_thresh_r = $y_thresh_l; $y_thresh_r <= $y_max;
                                            $y_thresh_r += $y_incr) {

         #__GET 2x2 TABLE CELL COUNTS
            my ($top_left, $bot_left, $top_right, $bot_right) = _counts_ (
               $pts, $x_thresh, $y_thresh_l, $y_thresh_r
            );

         #__PERFORM 2x2 TABLE SIGNIFICANCE TEST
            my ($pval, $grand_total) = _test_ ($top_left, $bot_left, $top_right, $bot_right);

         #__DIAGNOSTIC TALLYING
            unless ($eval_done) {
               $total_gene_evals += $grand_total;
               $eval_done = 1;
            }

         #__TERMINATE ONCE SIGNIFICANCE VANISHES: CURRENT Y-THRESHOLD IS...
            if ($pval > $pvalue_threshold) {

            #__OUTPUT
               my $p_val_typical_gene = 10**(-$y_thresh_l); #Log10
               my $p_val_large_gene = 10**(-$y_thresh_r);   #Log10
               print "   total gene evalluations for this cancer: $grand_total\n";
               print "d  GENE SIZE: TYPICAL GENES <= $x_thresh <= LARGE GENES";
            #  print "   GENE SIZE: TYPICAL GENES <= $x_thresh <= LARGE GENES";
               print "   (unchanged from default)" if
                     $x_thresh eq $x_thresh_gene_size;
               print "\n";
               print "   RECOMMENDED -LN(P-VALUE) CUT-OFFS\n";
               print "d     typical genes: -LN(P) = $y_thresh_l " .
            #  print "      typical genes: -LN(P) = $y_thresh_l " .
                     "(P-val = $p_val_typical_gene)";
               print "   (unchanged from default)" if
                     $y_thresh_l eq $y_thresh_fixed;
               print "\n";
               print "d     large genes:   -LN(P) = $y_thresh_r " .
            #  print "      large genes:   -LN(P) = $y_thresh_r " .
                     "(P-val = $p_val_large_gene)";
               print "   (unchanged from default)" if
                     $y_thresh_r eq $y_thresh_fixed;
               print "\n";
               print "   CALCULATION META-DATA\n";
               print "      total iterations = $iters\n";
               print "      final 2x2 table:   $top_left,  $top_right\n";
               print "                         $bot_left,  $bot_right\n";
               print "      final P-val $pval > threshold $pvalue_threshold " .
                     "(mutational significance status not statistically " .
                     "related to gene size)\n";

               my ($old_tps, $old_fps, $old_tns, $old_fns) = (0, 0, 0, 0);
               my ($new_tps, $new_fps, $new_tns, $new_fns) = (0, 0, 0, 0);
               my ($original_fp, $current_fp, $true_pos) = ({}, {}, {});
               foreach my $y (sort _numeric_ keys %{$gene_data}) {
                  foreach my $x (keys %{$gene_data->{$y}}) {
                     foreach my $gene_name (keys %{$gene_data->{$y}->{$x}}) {

                     #__ORIGINAL (FIXED-Y) FILTERING (NO LONG-GENE POST-HOC FILTER)
                        if ($y <= $y_thresh_l) {
                           if (exists $pos_genes->{$gene_name}) {
                              $old_fns++;
                           } else {
                              $old_tns++;
                           }
                        } else {
                           if (exists $pos_genes->{$gene_name}) {
                              $old_tps++;
                           } else {
                              $old_fps++;
                           }
                        }

                     #__NEW FILTERING RESULTS WITH LONG-GENE POST-HOC FILTER
                        if ($x > $x_thresh && $y > $y_thresh_l && $y < $y_thresh_r) {
                           $all_filtered_long_genes->{$gene_name}++;
                           $current_fp->{$gene_name} = "($x, $y)";
                           if (exists $pos_genes->{$gene_name}) {
                              $new_fns++;
                           } else {
                              $new_tns++;
                           }
                        } elsif ($y <= $y_thresh_l) {
                           $original_fp->{$gene_name} = "($x, $y)";
                           if (exists $pos_genes->{$gene_name}) {
                              $new_fns++;
                           } else {
                              $new_tns++;
                           }
                        } else {
                           $true_pos->{$gene_name} = "($x, $y)";
                           if (exists $pos_genes->{$gene_name}) {
                              $new_tps++;
                           } else {
                              $new_fps++;
                           }
                        }
                     }
                  }
               }
               print "      OLD FILTERING: fp = $old_fps fn = $old_fns\n";
               print "      NEW FILTERING: fp = $new_fps fn = $new_fns\n";
           #   print "      genes originally filtered as false positives (* = NOT part of 366 gene set)\n";
           #   foreach my $gene_name (sort keys %{$original_fp}) {
           #      print "d      ";
           #      if (exists $pos_genes->{$gene_name}) {
           #         print "  ";
           #      } else {
           #         print "* ";
           #      }
           #      print "$gene_name   $original_fp->{$gene_name}\n";
           #   }
               print "      additional genes now filtered as FPs under 'large gene' test (* = NOT part of 366 gene set)\n";
               foreach my $gene_name (sort keys %{$current_fp}) {
                  print "d      ";
                # print "       ";
                  if (exists $pos_genes->{$gene_name}) {
                     print "  ";
                  } else {
                     print "* ";
                  }
                  print "$gene_name   $current_fp->{$gene_name}\n";
               }
               print "      genes kept as significant (* = NOT part of 366 gene set)\n";
               foreach my $gene_name (sort keys %{$true_pos}) {
                  print "d      ";
                # print "       ";
                  if (exists $pos_genes->{$gene_name}) {
                     print "  ";
                  } else {
                     print "* ";
                  }
                  print "$gene_name   $true_pos->{$gene_name}\n";
               }
               print "\n";
               goto NEXT_FILE;
            } else {
               $iters++;
            }
         }
      }
      }
      print "\n     NO CONVERGENCE USING CURRENT PARAMETERS\n\n";
      NEXT_FILE:
   }

#__DIAGNOSTIC GRAND TALLY
   print "\n";
   print "TOTAL GENE EVALUATIONS OVER ALL $num_cancers CANCERS = $total_gene_evals\n\n";

#__UNION OF ALL GENES NETTED BY THE LONG GENE FILTER
   print "\nunionized list of all genes filtered as FPs under 'large gene' test over all cancers (* = part of 366 gene set)\n";
   foreach my $gene_name (sort keys %{$all_filtered_long_genes}) {
      if (exists $pos_genes->{$gene_name}) {
         print "* ";
      } else {
         print "  ";
      }
      print "$gene_name   ($all_filtered_long_genes->{$gene_name} cancers)\n";
   }

##################
#  POST PROCESS  #
##################

#   These are examples of "good" genes that should not be kicked out
#   gene     max -ln(p)  size
#   ARID2    7.67         5505
#   ATM      7.16         9168
#   KMT2A    7.07        11916
#   KMT2B    12.08        8145
#   NF1      6.775        8517
################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################

sub _numeric_ {$a <=> $b}

#   THE LINE "if (defined $y && $y) {" **IS** READING OCCURENCES OF 0.0
#   IN THE INPUT FILES BECAUSE, ALTHOUGH PERL PROCESSES NUMBERS AND
#   TEXT-THAT-RESEMBLES-A-NUMBER DIFFERENTLY (SEE BELOW EXAMPLE), IT *READS*
#   A FILE OF NUMBERS INITIALLY AS TEXT
#
#   bash-3.2$ perl             bash-3.2$ perl
#   $s = 0.0;                  $s = "0.0";   <---quotation marks on this one
#   if ($s) {                  if ($s) {
#     print "yes\n";              print "yes\n";
#   } else {                   } else {
#     print "no\n";               print "no\n";
#   }                          }
#   no                         yes

sub _read_file_ {
   my ($file) = @_;
   open (F, $file) || die "cant open $file";
   my ($pts, $gene_data, $max_y_large_gene) = ([], {}, 0);
   while (<F>) {

   #__PARSE LINE
      next if /^#/;
      chomp;
      my ($gene_name, $x, $y) = split;

   #__PROCESS A LEGITIMATE LINE
      if (defined $y && $y) {

      #__STORE POINTS AND TRACK MAX Y-VALUE FOR LARGE GENES AND GENE NAMES
         push @{$pts}, [$x, $y];
         if ($x > $x_thresh_gene_size) {
            $max_y_large_gene = $y if $y > $max_y_large_gene;
         }
         $gene_data->{$y}->{$x}->{$gene_name}++;
#        print "$x   $y\n";  # for xmgrace plotting
      }
   }
   close (F);
   return ($pts, $gene_data, $max_y_large_gene);
}

sub _counts_ {
   my ($pts, $x_thresh, $y_thresh_l, $y_thresh_r, $x_min) = @_;

#__DEFAULTS
   $y_thresh_r = $y_thresh_l unless $y_thresh_r;
   $x_min = 0 unless $x_min;
 
#__CREATE OBSERVED 2X2 TABLE THAT RESULTS FROM THESE THRESHOLD SETTINGS
   my ($top_left, $bot_left, $top_right, $bot_right) = (0, 0, 0, 0);
   foreach my $pt (@{$pts}) {
      my ($x, $y) = @{$pt};

   #__LEFT-HALF-PLANE FOR TYPICAL GENES COMPARES TO FIXED Y-THRESHOLD
      if ($x <= $x_thresh) {
         if ($y <= $y_thresh_l) {
            $bot_left++;
         } else {
            $top_left++;
         }

   #__AND RIGHT-HALF-PLANE FOR LARGE GENES TO PROVISIONAL Y-THRESHOLD
      } else {
         if ($y <= $y_thresh_r) {
            $bot_right++;
         } else {
            $top_right++;
         }
      }
   }
   return ($top_left, $bot_left, $top_right, $bot_right);
}

sub _test_ {
   my ($top_left, $bot_left, $top_right, $bot_right) = @_;

#__MARGINAL TOTALS FOR THESE THRESHOLDS
   my $row_top = $top_left + $top_right;
   my $row_bot = $bot_left + $bot_right;
   my $col_left = $top_left + $bot_left;
   my $col_righ = $top_right + $bot_right;

#__GRAND TOTAL
   my $grand_total = $bot_left + $top_left + $bot_right + $top_right;

#__SANITY CHECKING
#  die "tallying problem" unless $row_top + $row_bot == $grand_total;
#  die "tallying problem" unless $col_left + $col_righ == $grand_total;

#__THEORETICAL 2X2 TABLE OF EXPECTED VALUES FROM MARGINAL TOTALS
   my $top_left_expec = $col_left * $row_top / $grand_total;
   my $bot_left_expec = $col_left * $row_bot / $grand_total;
   my $top_right_expec = $col_righ * $row_top / $grand_total;
   my $bot_right_expec = $col_righ * $row_bot / $grand_total;

#__COMPUTE CHI-SQUARE STATISTIC
   my $chisq = ($top_left - $top_left_expec)**2 / $top_left_expec +
               ($top_right - $top_right_expec)**2 / $top_right_expec +
               ($bot_left - $bot_left_expec)**2 / $bot_left_expec +
               ($bot_right - $bot_right_expec)**2 / $bot_right_expec;

#__2X2 TABLE TEST HAS DOF OF 1
   my $dof = 1;

#__COMPUTE GOODNESS-OF-FIT P-VALUE OF OBSERVED VS EXPECTED
   my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);
   return ($pval, $grand_total);
}

################################################################################
#                                                                              #
#                                   D A T A                                    #
#                                                                              #
################################################################################

#  EXPLANATION OF DATA: LETS US ESTIMATE SENSITIVITY & SPECIFICITY
#
#  THESE ARE THE GENES WHICH SHOULD BE KEPT IN THE ANALYSIS --- IT IS AN ERROR
#  TO FILTER THEM --- ALL OTHER GENES SHOULD BE FILTERED

__DATA__
TP53
PIK3CA
BRAF
PTEN
LRP1B
RYR2
APC
KRAS
SYNE1
ARID1A
ZFHX4
DNAH5
IDH1
PCLO
KMT2C
HMCN1
SPTA1
FAT3
FAT4
KMT2D
RYR3
EGFR
CTNNB1
USH2A
ATRX
NF1
DMD
FAT1
VHL
RELN
DNAH9
COL11A1
GPR98
LRP2
PBRM1
NAV3
RB1
ATM
NOTCH1
FBXW7
CDKN2A
CSMD2
DNAH11
PTPRD
HUWE1
CACNA1E
FBN2
NRAS
TENM1
PIK3R1
BAI3
DNAH7
NALCN
NRXN1
SETD2
PTPRT
ANK2
RIMS2
SACS
ANK3
LRRC7
CREBBP
SMARCA4
CDH1
CNTNAP2
DNAH2
CTNNA2
DOCK2
TAF1L
PCDH10
ZFHX3
CHEK2
EP300
PREX2
SMAD4
ERBB4
KMT2A
MED12
SVEP1
KEAP1
BRINP3
EPHA5
ZNF536
CHD4
ARID2
NBEA
SCN1A
GRIN2A
APOB
CACNA1C
ASTN1
SCN3A
NIPBL
CDH9
THSD7B
CHD7
CTNND2
LRFN5
GATA3
LPHN3
LRRK2
NRXN3
PLXNA4
UNC79
FREM2
KDM6A
MTOR
MYH4
PCDH17
TENM2
ROBO2
TENM3
HCN1
KCNT2
EPHA3
SRCAP
CDH8
HECW1
PCDH11X
SPEN
TENM4
ABCA12
DNMT3A
TRPS1
FRAS1
TSHZ3
ADCY2
DPP10
MAP3K1
NCOR1
SLIT2
SORCS1
CDH18
NBPF1
NSD1
MYH9
PCDH9
ADAMTS20
NCOR2
USP9X
CIC
KMT2B
NFE2L2
ZNF521
ASXL3
DSP
STK11
COL22A1
HYDIN
ATR
CPS1
DOCK4
SCN5A
SCN9A
TAF1
COL1A2
NEFH
STAG2
BCOR
SPTAN1
ARHGAP35
CACNA1D
CDH12
PLCB1
ABCC9
ADAMTS12
CDC27
EPHB1
GRM1
ZFPM2
AFF2
EP400
SMG1
BAP1
COL3A1
CTCF
GRM3
OTOF
SALL1
ASTN2
CACNA1A
EPHA6
THSD7A
BSN
ADAMTS16
BCLAF1
GRIA2
GRID1
PAPPA
WBSCR17
ABCB1
CHD8
COL5A2
SF3B1
BAZ2B
COL4A5
DOCK10
FLT3
FSTL5
LPHN2
NOTCH3
SLITRK5
TRHDE
ZEB2
ADAMTS19
ERBB3
GRM8
LRP6
CNKSR2
DOCK3
SLITRK2
ALK
ERBB2
PPFIA2
FMN2
GRIK2
GRM5
PLCL1
DCLK1
KDR
KIAA2022
KCNN3
KDM5C
PLCB4
SLC8A1
ARHGAP5
PCDH18
SEMA5A
ANO4
CNTN5
DLC1
PTPRB
EPB41L3
LRRC4C
MGAM
SETBP1
THOC2
ZNF292
CHD9
ATRNL1
PABPC3
RIMBP2
TUBA3C
FRYL
NDST4
PCDH19
MECOM
YLPM1
BCORL1
CNTN6
HDAC9
KCNH8
KLHL4
NLGN4X
PHF2
SLC12A5
TLL1
CACNA2D1
DOCK8
ITGA8
OTOGL
SPOP
TIAM1
TRPC4
CHD5
GRIA3
GRIN3A
KAT6B
KCNH1
LTBP1
MMP16
RUNX1T1
ST18
BRINP1
SIPA1L2
CDK12
CHD3
EPHA7
HGF
NTRK3
PDGFRA
SCN11A
SEMA6D
SLC4A4
TSHZ2
CARD11
CDH6
NCAM2
PPP2R1A
TNIK
TNN
ZNF804A
AFF3
CASZ1
GPR158
KCND2
KCNQ3
KSR2
RBFOX1
ZMYM3
ZNF423
ATP11C
FOXP2
LINGO2
MAP3K4
MET
PDZRN4
TNRC18
FGFR2
PCDH7
PTCHD4
RASA1
RERE
SMC1A
SORL1
SOX5
ADAMTS2
ARFGEF2
BZRAP1
CASP8
CLIP1
FBXL7
GALNT13
NLGN1
PDE1C
SLITRK4
ANKRD11
CDH11
DPYD
ELMO1
GIGYF2
HRAS
RUNX1
ASXL1
CADPS2
GABRA6
NF2
PXDNL
ACTN2
ELTD1
RBM10
TP63
UBBP4
ATXN1
GABRA1
INPPL1
JAK2
SLIT3
SMARCA1
DOCK5
HEPH
KCNQ5
MAP1B
NPM1
RGPD4
TBX3
TSC2
COL5A3
DICER1
LPPR4
PPP1R9A
PTCH1
SEMA3A
ACAN
ANKRD50
BCL11A
FRG1
KCNJ12
PCSK5
SALL3
SATB2
WSCD2
LZTR1
SCAND3
DGKI
DYNC1I1
NETO1
PLCG2
ZEB1
ATP1A2
MAP2K4
PKP4
SHANK1
TCF7L2
XYLT1
ZNF384
ARID5B
C6
CUL3
GPC6
MAGI1
MN1
RPS6KA6
SFMBT2
SGIP1
TET2
UNC5C
ZFP36L2
ZNF711
AMER1
CNTN4
FAM184A
FOLH1
FOXA1
H3F3A
INO80
LAMB2
PDZRN3
POLR2A
SIN3A
SLITRK6
SMC3
SOX9
VAV3
WT1
BCHE
LURAP1L
AR
ASXL2
CUL1
ENPP2
GABRB3
KANSL1
KIAA1211
MAP3K9
PDS5B
TGFBR2
EPHA2
GUCY1A2
MCC
MED1
EEF1A1
PARG
ZBTB20
BRCA1
DDX3X
KIAA0907
RANBP6
SEMA3D
SKIDA1
TRPV6
ZNF429
ARHGEF6
HSPA8
NTM
NTNG1
OGT
PCF11
PTPN11
RBM47
ADAM23
CNTNAP1
DENND5B
KLHL13
MED13L
MED15
MGAT4C
TBP
USP7
CD163
COL21A1
GCC2
ITGA4
KCNA1
TRERF1
TRPC3
ZNF385D
ABL1
ANO2
CHD1
ESRRG
JAK1
KDM6B
NEBL
PAK7
SDHA
AKT1
COL4A4
IL32
INHBA
RGPD3
ACVR1
FUBP1
RBM5
ACVR2A
ATXN2
AXIN1
DAAM1
ELAVL2
HSP90AB1
LATS1
LPAR4
MST1
MYO1B
TRANK1
VPS8
ZFP36L1
BMP2K
CSDE1
FAM155A
FAP
HIST1H3B
IKZF1
MB21D2
PAK3
RNF43
SEMA3C
DACH1
DKK2
DNM2
DTX1
ELFN2
FGFR3
GLUD2
HDX
JAG1
MAML2
PLA2G4A
PODXL
RBM26
TCF12
TMTC1
TVP23C-CDRT4
EZH2
HAS2
IDH2
MYB
PAX5
PLEKHA5
POTEE
SNCAIP
BICD2
BRD4
DIS3
DYRK1A
FAM91A1
NAA15
PIK3C2G
PLCG1
RAD21
TMCC1
TUBB8
XPO1
AJUBA
DBR1
INPP4A
LDB2
NOS3
PLXDC2
RBMXL1
STRN
SUPT16H
WAC
CBL
CEP135
CRIPAK
HS6ST1
IRF2BPL
KCNA5
KCTD16
NUP107
PGR
POLG
RHOA
TSC1
AMY2B
DDX42
FOXP1
KCNG1
MEF2A
PABPC5
RPL5
SP8
ZNF800
DDX5
JAK3
KCNK2
KIT
LIMCH1
LIMK2
SLC9A5
SMARCB1
TOX
U2AF1
USP8
ELF3
ERCC2
ETV6
FZD10
KCNA3
MPP7
NEGR1
PCMTD1
RPS6KA3
SLC6A13
STK19
TLE1
ACTG1
BHLHE22
CALN1
CDKN1B
CLVS1
FIP1L1
GPX1
KLF5
LCP1
MEN1
NEFM
SCAF11
UBC
USP51
VGLL3
WDR47
ACTB
ASB2
CSGALNACT2
HNF1A
KBTBD6
LINGO1
NOVA1
NUP93
PCSK1
RAC1
RXRA
TXNIP
ZNF717
ZNRF3
CNTNAP3B
EOMES
ESR1
FXR2
GALNT2
HNRNPM
MAP2K3
SLC5A3
SUZ12
ZC3H11A
ALB
BCL6B
FAM46A
HSPA1L
MORC4
PHF6
RARG
VAV2
ALCAM
CTBS
EFEMP1
GHR
HNRNPL
HTR7
MGAT3
MYC
PHACTR1
PPM1D
RPL22
TUSC3
ZBTB33
CBFB
CLASRP
FOXA2
GRIN2D
GUF1
KCNJ4
LRRIQ3
LUM
MAMLD1
MAPK1
PDE1A
PLOD2
SLC1A3
SLC25A5
SMAP1
ARVCF
B2M
IL7R
KAT2A
LRRC4
METTL14
OLFML2A
RALY
TMPRSS6
ATP6V1C2
DOK6
EIF1AX
EIF2B5
FAM65A
FEM1A
HOXD10
PBX2
PIK3R5
SPRED1
SRRM1
TP73
WIF1
ABCD1
CHST1
DDX39B
FOSL2
GATA2
GNAQ
GPS2
KRTAP5-5
MAX
MBD1
MTA2
PPP6C
TAB2
ACADS
CDKN1A
ELF1
G6PD
IL6ST
KDM1A
KDM4C
KRT5
MYCN
PAIP1
PAK2
PCBP1
PHLDA1
ZBTB7B
CAMK1D
HDAC2
HSPA6
LRRC55
MAP7
QKI
SGK1
SLC20A1
SUFU
ZC3H15
ZRANB2
ZSWIM6
AGPS
AGTR1
CALD1
CEL
CLCN7
CNOT11
GNB1
HIST1H2BD
HNRNPA1
NUMB
ONECUT2
P2RY1
PRKAR1A
SOX17
VMP1
AGGF1
ANKRD34A
DCAF13
HSPBP1
LMX1B
MUC12
NMT2
PCDHB16
PSPH
PTRF
RARA
RFX5
SCAF1
ZNF513
CCDC41
CDK11B
FZD9
GRAMD1A
KIF7
LNX2
NCK1
NPNT
NPY1R
TDG
TOB1
WDR45
ZBED4
ZNF516
ZXDB
CDK11A
GAL3ST3
HNRNPK
IK
NET1
PEX5
SAV1
SIAH1
SMTNL2
ZNF43
ZNF750
AGAP6
CEBPA
CNPY3
FAM120B
FZD8
MAP3K3
MSI1
NIPA2
PEX14
USP10
ZBTB39
ZIC5
ACADM
AOC3
CALR
DKFZP761J1410
ETV3
FAM117B
FOXD4
GATA6
NACA2
PDIA4
RBBP4
RP11-14N7.2
WASL
ZNF706
ADRA2B
ANKLE1
BCL2L11
CCND1
CERS2
FBXW11
GNAL
HRCT1
RBFOX2
RRN3
SIRPA
SOX4
TGIF1
ATP1B1
B3GNT6
C1GALT1C1
DYNC1I2
FASTKD5
GGT1
PMFBP1
SH2B3
SNRNP70
ZNF41
AZIN1
EPHX4
ERCC6L
HOXB13
HOXD8
HPD
HPS6
LRFN4
MARCKS
NKX2-1
PTP4A3
PUS1
RRAD
RRS1
SUCLG1
ARFGAP3
CDKN2C
FAM57B
JUNB
PAQR4
PYGO2
RHOB
SEPT7
SLA
SLC47A2
TAL1
CCNB2
FBXO46
GATA1
HIST1H3D
KRT8
PKDREJ
PPIAL4G
RRAS2
ATG5
CHMP1B
FOXQ1
GCNT3
IGHV3OR16-8
NKPD1
PRSS27
RNPC3
RPL14
RPP25L
RYBP
SSTR5
ZNF469
CTDNEP1
EXOSC2
FAM186A
GOLGA6L19
HOXC8
MXD4
NKX3-1
PRSS22
SLC25A45
ABHD17A
F12
GATSL3
HCFC1R1
HTRA2
KLK1
KRTAP5-1
NHLRC1
RIMBP3
SH3BP2
SLMO2
TSPAN2
VEGFB
ZNF697
CD7
CREG2
GADD45GIP1
MTHFD2L
SERF2
SULT1C4
TCEB1
TOB2
ADRA2A
AL592183.1
CDK3
FUT6
GSG2
HEBP1
INPP5E
NAT9
RBKS
CEBPB
LYPLA2
RAB27A
ABHD14A
ADO
CBX7
FAM104B
IL8
MESDC1
NRARP
RASSF7
AC026703.1
AC090616.2
AL133373.1
ASIP
NDUFA11
SLPI
