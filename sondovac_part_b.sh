#!/bin/bash

# Determine script's directory
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load functions shared by both parts, introductory message
source $SCRIPTDIR/sondovac_functions || {
  echo
  echo "Fatal error! Unable to load file \"sondovac_functions\" with required functions!"
  echo "It must be in same directory as \"$0\""
  echo "Check it and if needed download again whole script from"
  echo "https://github.com/V-Z/sondovac/"
  echo
  exit 1
  }

echo "This is part B of the pipeline."
echo
echo "This part processes assembly output from Geneious and produces the final list of low-copy nuclear probe sequences."

# Default values
# Counter if not both -i and -n options are used
CHECKMODE=0
# If not specifying explicitly otherwise (using -n), running in interactive mode
STARTINI="I"

# Parse initial arguments
while getopts "hvulrpeinf:b:d:" START; do
  case "$START" in
    h|v)
      echo "Usage options:"
      echo -e "\t-h, -v\tPrint this help message and exit"
      echo -e "\t-u\tCheck for updates and download them if needed and confirmed by user, and exit"
      echo -e "\t-l\tDisplay LICENSE for license information (this script is licensed under GNU GPL v.3, other software under variable licenses) and exit"
      echo -e "\t-r\tDisplay README for detailed usage instructions and exit"
      echo -e "\t-p\tDisplay INSTALL for detailed installation instructions and exit"
      echo -e "\t-e\tDisplay detailed citation information and exit"
      echo -e "\t-i\tRunning in interactive mode - script will on-demand ask for required input files, installation of missing software etc."
      echo -e "\t\tThis is recommended default value (the script runs interactively without explicit using option ${BOLD}-n${NORM})."
      echo -e "\t-n\tRunning in non-interactive mode. User ${BOLD}must${NORM} provide at least five input files below:"
      echo -e "\tYou can use ${BOLD}only one${NORM} of parameters ${BOLD}-i${NORM} or ${BOLD}-n${NORM} (not both of them)"
      echo
      echo -e "\tIf option ${BOLD}-f${NORM} is used and script is running in interactive mode, this value will be used as defaults, but may be later overridden."
      echo
      echo -e "\tOptions required for running in non-interactive mode:"
      echo -e "\t-f\tInput file in FASTA format (output of Geneious assembly)"
      echo
      echo -e "\tOther optional arguments (if not provided, default values are used):"
      echo -e "\t-b\tBait length"
      echo -e "\t\tDefault value: 120 (optimal length for phylogeny, use integer between 120 and 200)"
      echo -e "\t-d\tSequence similarity between the developed probe sequences (parameter \"-c\" of cd-hit-est, see its manual for details)"
      echo -e "\t\tDefault value: 0.9 (use decimal number ranging from 0 to 1)"
      echo -e "\t${BOLD}WARNING!${NORM} If parameters ${BOLD}-b${NORM} or ${BOLD}-d${NORM} are not provided, default values are taken and it is not possible to change them later (not even in interactive mode)."
      echo
      exit 2
      ;;
    u)
      scriptupdater
      ;;
    l)
      licenser
      ;;
    r)
      readmeview
      ;;
    p)
      installview
      ;;
    e)
      citationreference
      ;;
    i)
      echo "Running in interactive mode..."
      STARTINI="I"
      CHECKMODE=$((CHECKMODE+1))
      ;;
    n)
      echo "Running in non-interactive mode..."
      STARTINI="N"
      CHECKMODE=$((CHECKMODE+1))
      ;;
    f)
      SEQUENCES=$OPTARG
      echo "Input file: $SEQUENCES"
      ;;
    b)
      BAITL=$OPTARG
      if [[ "$BAITL" =~ ^[0-9]+$ ]] && [ "$BAITL" -ge 120 -a "$BAITL" -le 200 ]; then
	echo "Bait length: $BAITL"
	else
	  echo
	  echo "Error! For parameter \"-b\" you did not provide an integer of range 120-200!"
	  echo
	  exit 1
	fi
      ;;
    d)
      CDHITSIM=$OPTARG
      if [ "$(echo 0 '<=' $CDHITSIM | bc -l)" = 1 ] && [ "$(echo $CDHITSIM '<=' 1 | bc -l)" = 1 ]; then
	echo "Sequence similarity: $CDHITSIM"
      else
	echo
	echo "Error! For parameter \"-d\" you did not provide decimal number ranging from 0 to 1!"
	echo
	exit 1
      fi
      ;;
    ?)
      echo "Invalid option(s)!"
      echo "See \"$0 -h\" for usage options."
      exit 1
      ;;
    esac
  done

if [ "$CHECKMODE" == 2 ]; then
  echo
  echo "${BOLD}Error!${NORM} You provided both parameters ${BOLD}-i${NORM} (interactive mode) and ${BOLD}-n${NORM}"
  echo "(non-interactive mode). You ${BOLD}may${NORM} use ${BOLD}only one${NORM} of them (${BOLD}either -i${NORM} or ${BOLD}-n${NORM})!"
  echo
  exit 1
  fi

# Check which operating system the script is running on

# Ensure user reads introductory information
confirmgo

# Warn user this is not yet for daily usage
devrelease

# Check operating system
oscheck

# Set variables for working directory and PATH
workdirpath

# Check availability of all needed binaries

# Check if grep is available
checktools grep

# Check if egrep is available
checktools egrep

# Check if geneious_column_separator.pl is available
# checktools geneious_column_separator.pl

# Check if awk is available
checktools awk

# Check if sed is available
checktools sed

# Check if join is available
checktools join

# Check if cat is available
checktools cat

# Function to compile CD-HIT
function compilecdhit {
  {
  checktools make &&
  checktools g++ &&
  echo &&
  echo "Compiling CD-HIT from source code" &&
  echo &&
  cd $1 &&
  make -s openmp=yes || make -s &&
  cp -a *.pl $BIN/ &&
  cp -a cd-hit* $BIN/ &&
  cd $WORKDIR &&
  echo &&
  echo "\"CD-HIT\" is available. OK."
  } || {
    echo
    echo "Compilation failed. Please go to https://github.com/weizhongli/cdhit,"
    echo "download cd-hit-*.tgz, compile it and ensure it is in PATH."
    echo "Check last error messages to find why compilation failed."
    echo
    exit 1
    }
  }

# Check if cd-hit-est is available
{ command -v cd-hit-est >/dev/null 2>&1 && echo "\"cd-hit-est\" is available. OK."; } || {
  echo "\"cd-hit-est\" is required but not installed or available in PATH"
  echo
  echo "Type \"C\" to compile CD-HIT 4.6.1 from source available together with this script"
  echo "Type \"S\" to download latest CD-HIT source from https://github.com/weizhongli/cdhit and compile it"
  echo "Type \"B\" to copy CD-HIT 4.6.1 binary available together with the script (available for Linux and Mac OS X)"
  echo "Type \"M\" for manual installation - script will exit and you will have to install CD-HIT yourselves"
  read CDHIT
  while :
  do
    case "$CDHIT" in
      C|c)
	compilecdhit $SCRIPTDIR/src/cd-hit-v4.6.1-2012-08-27
	break
	;;
      S|s)
	checktools curl
	checktools unzip
	curl -L -o cd-hit-master.zip https://github.com/weizhongli/cdhit/archive/master.zip
	unzip -nq cd-hit-master.zip
	compilecdhit cdhit-master
	break
	;;
      B|b)
	echo "Copying CD-HIT binaries"
	case "$OS" in
	  Mac)
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/cd-hit* $BIN/
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/clstr*.pl $BIN/
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/plot_*.pl $BIN/
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/psi-cd-hit*.pl $BIN/
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/make_multi_seq.pl $BIN/
	    ;;
	  Linux)
	    if [ "$OSB" == "64b" ]; then
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/cd-hit* $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/clstr*.pl $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/plot_*.pl $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/psi-cd-hit*.pl $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/make_multi_seq.pl $BIN/
	    else
	      cp -p $SCRIPTDIR/pkgs/linux32b/bin/cd-hit* $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux32b/bin/clstr*.pl $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux32b/bin/plot_*.pl $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux32b/bin/psi-cd-hit*.pl $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux32b/bin/make_multi_seq.pl $BIN/
	    fi
	    ;;
	  *) echo
	    echo "Binary is not available for $OS $OSB. Going to compile it from source code."
	    echo
	    compilecdhit $SCRIPTDIR/src/cd-hit-v4.6.1-2012-08-27
	    ;;
	esac
	break
	;;
      M|m) echo && echo "Please, go to http://weizhongli-lab.org/cd-hit/ and install CD-HIT and ensure it is in PATH" && echo && exit;;
      *) echo "Wrong option. Use C, S, B or M." && read CDHIT;;
    esac
  done
  }
echo "Rest not finished yet" && exit
# Variables

# Geneious output files are infiles here - statistics of contigs
TSVLIST=$1
# Geneious output files are infiles here - consensus and unused sequences
SEQUENCES=$2
# Sequences converted from FASTA to tabular format
SEQUENCESTAB="${SEQUENCES%.*}.tsv"


# Part 3: Assemble the obtained sequences in contigs (part B)

# Check if TSV output of Geneious contains at least requested columns
echo
if grep -q "# Sequences" $TSVLIST; then
    echo "Column \"# Sequences\" is presented in $TSVLIST. OK."
    if grep -q "Pairwise Identity" $TSVLIST; then
      echo "Column \"Pairwise Identity\" is presented in $TSVLIST. OK."
      if grep -q "Description" $TSVLIST; then
	echo "Column \"Description\" is presented in $TSVLIST. OK."
	if grep -q "Mean Coverage" $TSVLIST; then
	  echo "Column \"Mean Coverage\" is presented in $TSVLIST. OK."
	  if grep -q "Name" $TSVLIST; then
	    echo "Column \"Name\" is presented in $TSVLIST. OK."
	    if grep -q "Sequence Length" $TSVLIST; then
	      echo "Column \"Sequence Length\" is presented in $TSVLIST. OK."
	    else
	      echo "Error! Column \"Sequence Length\" is missing!"
	      echo "Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\". Please, export the TSV file again."
	      exit 1
	    fi
	  else
	    echo "Error! Column \"Name\" is missing!"
	    echo "Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\". Please, export the TSV file again."
	    exit 1
	  fi
	else
	  echo "Error! Column \"Mean Coverage\" is missing!"
	  echo "Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\". Please, export the TSV file again."
	  exit 1
	fi
      else
	echo "Error! Column \"Description\" is missing!"
	echo "Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\". Please, export the TSV file again."
	exit 1
      fi
    else
      echo "Error! Column \"Pairwise Identity\" is missing!"
      echo "Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\". Please, export the TSV file again."
      exit 1
    fi
  else
    echo "Error! Column \"# Sequences\" is missing!"
    echo "Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\". Please, export the TSV file again."
    exit 1
  fi
echo

# Check if TSV output of Geneious contains only required columns or more
if egrep -q "# Sequences[[:blank:]]+% Pairwise Identity[[:blank:]]+Description[[:blank:]]+Mean Coverage[[:blank:]]+Name[[:blank:]]+Sequence Length" $TSVLIST
  then
    echo "OK, $TSVLIST is correct input file."
    TSVLIST2=$TSVLIST
  else
    echo "Input file $TSVLIST seems to contain more columns than required. Needed columns will be extracted."
    geneious_column_separator.pl $TSVLIST || echo "Extraction failed. Please, do it manually. Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\"."
    TSVLIST2="${TSVLIST%.*}.columns.tsv"
    echo "File with extracted columns was saved as $TSVLIST2 for possible later usage"
  fi
echo

# Check the statistics
# Check total number of bp
echo "Total number of base pairs:"
{ cut -f6 $TSVLIST2 | awk '$1>119' | awk '{s+=$1}END{print s}'; } || { echo && echo "${BOLD}Error!${NORM} XXX failed. Aborting." && echo && exit 1; }
# Check number of contigs
echo "Number of contigs:"
{ cut -f6 $TSVLIST2 | awk '$1>119' | wc -l; } || { echo && echo "${BOLD}Error!${NORM} XXX failed. Aborting." && echo && exit 1; }
echo

# Convert FASTA to TSV
echo "Converting FASTA to TAB"
# REWRITE!!!
perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' $SEQUENCES > $SEQUENCESTAB || { echo && echo "${BOLD}Error!${NORM} XXX failed. Aborting." && echo && exit 1; }
echo "Converted tabular format was saved as $SEQUENCESTAB for possible later usage" || { echo && echo "${BOLD}Error!${NORM} XXX failed. Aborting." && echo && exit 1; }


# Separate the assembled and unassembled sequences:
grep 'oxalis' oxalis.all-lengths.assembled+unassembled.lessthan1kbtranscript-hits_FINAL.tab > oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab
# 41952 contigs (Aaron: 41952)
grep -v 'oxalis' oxalis.all-lengths.assembled+unassembled.lessthan1kbtranscript-hits_FINAL.tab > oxalis.all-lengths.unassembled.lessthan1kbtranscript-hits_FINAL.tab
# 46061 unassembled sequences (Aaron: 46061)
# Filter the file with the assembled sequences – count the assemblies (the ones indicated with "Contig") making up genes of ≥960 bp / ≥600 bp, comprised of putative exons ≥120 bp
awk '{print $1"\t"length($2)}' oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>959' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>959' | wc -l
# 473 genes of ≥960 bp (exons ≥120 bp), total 638816 bp (Aaron: 473, total 638816 bp)
awk '{print $1"\t"length($2)}' oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
# 1223 genes of ≥600 bp (exons ≥120 bp), total 1195988 bp (Aaron: 1223, total 1195988 bp)
awk '{print $1"\t"length($2)}' oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' > oxalis.probes.120+600bp_FINAL.txt
# Filter the file with the unassembled sequences:
awk '{print $1"\t"length($2)}' oxalis.all-lengths.unassembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>399' | wc -l
# 33 unassembled sequences making up genes of ≥400 bp (Aaron: 32)
awk '{print $1"\t"length($2)}' oxalis.all-lengths.unassembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
awk '{print $1"\t"length($2)}' oxalis.all-lengths.unassembled.lessthan1kbtranscript-hits_FINAL.tab | sed s/_/\\t/g | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
# 1 unassembled sequence making up a gene of ≥600 bp (exons ≥120 bp), total 669 bp (Aaron: 1, total 669 bp)

# Part 4: Create the final .fasta file for the Hyb-Seq probes.

# Extract and sort the assemblies making up genes of ≥600 bp:
sed s/^/Assembly_/ oxalis.probes.120+600bp_FINAL.txt | cut -f1 -d' ' | sort -k1,1 > oxalis.probes.120+600bp_FINAL.forjoin
# Make a file with all exons ≥120 bp (as the sequence labels of my GENEIOUS file oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_FINAL.tab looks slightly different from Aaron's, I changed them to his notation by replacing "no_N-1kb_hits_FINAL" with "no_N.no_1kb-hits" in a text editor and named the file oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_AaronsFormat_FINAL.tab):
awk '{print $1"\t"length($2)"\t"$2}' oxalis.all-lengths.assembled.lessthan1kbtranscript-hits_AaronsFormat_FINAL.tab | awk '$2>119' > oxalis.120bp+.assembled.lessthan1kbtranscript-hits_AaronsFormat_FINAL.tab
# Make the assembly number the first field and sort (I modified Aaron's command as it didn't work for me):
sed s/oxalis.blat.no_N.no_1kb-hits_//oxalis.120bp+.assembled.lessthan1kbtranscript-hits_AaronsFormat_FINAL.tab | sed s/_C/\\tC/ | sort -k1,1 > oxalis.120bp+.assembled.lessthan1kbtranscript-hits_AaronsFormat_sorted_FINAL.tab
# Make a file with all exons ≥120 bp and all assemblies making up genes of ≥600 bp:
join oxalis.probes.120+600bp_FINAL.forjoin oxalis.120bp+.assembled.lessthan1kbtranscript-hits_AaronsFormat_sorted_FINAL.tab > oxalis.probes.120+600bp_FINAL.tab
# 5232 assemblies making up genes of ≥600 bp (exons ≥120 bp) (Aaron: 5232)
# Convert .tab to .fasta (I modified Aaron's set of commands as it didn't work for me):
sed s/' '/_/ oxalis.probes.120+600bp_FINAL.tab | sed s/' '/_/ > oxalis.probes.120+600bp_modified_FINAL.tab
sed s/^/'>'/ oxalis.probes.120+600bp_modified_FINAL.tab | sed s/' '/\\n/ > oxalis.probes.120+600bp.assembled_FINAL.fsa
# IMPORTANT: When I created the file of assemblies making up genes of ≥600 bp oxalis.probes.120+600bp_FINAL.forjoin, only the assemblies indicated by "Contig" were listed. Therefore, the remaining assemblies from oxalis.120bp+.assembled.lessthan1kbtranscript-hits_Aarons-Format_FINAL.tab have to be selected and added to the .fasta file of the probes:
grep -v Contig oxalis.120bp+.assembled.lessthan1kbtranscript-hits_AaronsFormat_FINAL.tab | awk '$2>599' | sed s/^/'>'/ | sed s/\\t/_/ | sed s/\\t/\\n/ > oxalis.probes.120+600bp.contig0_FINAL.fsa
# 16 assemblies (Aaron: 16)
# Combine the two .fasta files:
cat oxalis.probes.120+600bp.assembled_FINAL.fsa oxalis.probes.120+600bp.contig0_FINAL.fsa > oxalis.probes.120+600bp_FINAL.fasta
# 5248 assemblies making up genes of ≥600 bp (exons ≥120 bp), total 1207375 bp (Aaron: 5248 assemblies, 1207375 bp)
# Note that for further file processing you need to modify the label of the 16 assemblies not indicated by "Contig". Remove "oxalis.blat.no_N.no_1kb-hits" and add "Contig_0" at the right position (between assembly number and assembly length). I did this in a text editor and named the file oxalis.probes.120+600bp.fasta.

# Part 5: Make the final quality control of the probe sequences before sending them to company for bait synthesis

# Check for sequence similarity between the developed probe sequences with CD-HIT-EST:
cd-hit-est -i oxalis.probes.120+600bp.fasta -o oxalis.probes.120+600bp.SimilarityTest_FINAL -c 0.9
# One of the three outfiles is a .fasta file: oxalis.probes.120+600bp.SimilarityTest_FINAL
# Adjust the number of assemblies making up genes of ≥600 bp as numerous assemblies were eliminated in the sequence similarity check. Consequently, there will be assemblies in the file oxalis.probes.120+600bp.SimilarityTest_FINAL which will not sum up to genes of ≥600 bp.
# Rename oxalis.probes.120+600bp.SimilarityTest_FINAL to oxalis.probes.120+600bp.SimilarityTest_FINAL.fasta.
# This file has to be converted to .txt:
fasta2tab.pl oxalis.probes.120+600bp.SimilarityTest_FINAL.fasta > oxalis.probes.120+600bp.SimilarityTest_FINAL.txt
# Count all assemblies, comprised of putative exons ≥120 bp:
awk '{print $1"\t"length($2)}' oxalis.probes.120+600bp.SimilarityTest_FINAL.txt | awk '{s+=$2;c++}END{print s}'
# total 1166757 bp of assemblies (also including the ones making up genes of <600 bp)
# Count the assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp:
awk '{print $1"\t"length($2)}' oxalis.probes.120+600bp.SimilarityTest_FINAL.txt | sed s/_/\\t/g | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' oxalis.probes.120+600bp.SimilarityTest_FINAL.txt | sed s/_/\\t/g | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
# 1168 genes of ≥600 bp (exons ≥120 bp), total 1131805 bp
awk '{print $1"\t"length($2)}' oxalis.probes.120+600bp.SimilarityTest_FINAL.txt | sed s/_/\\t/g | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' > oxalis.probes.120+600bp_2ndFINAL.txt
# Extract and sort the assemblies making up genes of ≥600 bp:
sed s/^/Assembly_/ oxalis.probes.120+600bp_2ndFINAL.txt | cut -f1 -d' ' | sort -k1,1 > oxalis.probes.120+600bp_2ndFINAL.forjoin
# Modify the assembly number and sort:
sed s/_C/\\tC/ oxalis.probes.120+600bp.SimilarityTest_FINAL.txt | sort -k1,1 > oxalis.probes.120+600bp.SimilarityTest_sorted_FINAL.tab
# Make a file with all exons ≥120 bp and all assemblies making up genes of ≥600 bp:
join oxalis.probes.120+600bp_2ndFINAL.forjoin oxalis.probes.120+600bp.SimilarityTest_sorted_FINAL.tab > oxalis.probes.120+600bp_2ndFINAL.tab
# Convert .tab to .fasta:
sed s/' C'/_/ oxalis.probes.120+600bp_2ndFINAL.tab | sed s/ontig/Contig/ | sed s/^/'>'/ | sed s/' '/\\n/ > oxalis_HybSeqProbes.fasta
# 4930 assemblies making up genes of ≥600 bp (exons ≥120 bp), total 1131805 bp
# This is how I calculated the total number of base pairs:
fasta2tab.pl oxalis_HybSeqProbes.fasta > oxalis_HybSeqProbes.txt
awk '{print $1"\t"length($2)}' oxalis_HybSeqProbes.txt | awk '{s+=$2;c++}END{print s}'
# Removing remaining cp/mt genes from probe set
blat -t=dna -q=dna -out=pslx Ricinus_communis_reference.fsa oxalis_HybSeqProbes.fasta blat_targetEnrichment_probeSequences_plastidSequences.pslx
# The resulting transcripts, such as below, should then be removed from the probe set.
# In a last step repetitive regions have to be masked. I used REPEATMASKER (http://www.repeat-masker.org/cgi-bin/WEBRepeatMasker).
# For a detailed introduction to REPEATMASKER see http://sebastien.tempel.free.fr/Boulot/UsingRepeat-Masker.pdf. I used the following settings:
# Search engine: cross_match (often more sensitive than the other engines)
# Speed / sensitivity: slow
# DNA source: other, please specify; Populus trichocarpa
# Return format: html
# Return method: email
# Masking options: repetitive sequences replaced by strings of N
# Repeat options: mask interspersed and simple repeats
# Matrix: RepeatMasker choice
# Artifact check: skip bacterial insertion element check
# Amongst the four results files there is a file with masked repetitive regions in .fasta (so-called masked file). Send this file (I renamed it to oxalis_HybSeqProbes_repeatsmasked.fasta) to MYcroarray for probe synthesis.

exit
