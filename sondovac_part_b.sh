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
# Bait length
BAITL=120
# CD-HIT sequence similarity
CDHITSIM=0.9
# Create empty variables for file names
REFERENCECP=""
TSVLIST=""
SEQUENCES=""

# Parse initial arguments
while getopts "hvulrpeinc:x:z:b:d:" START; do
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
      echo -e "\tIf options ${BOLD}-c${NORM}, ${BOLD}-x${NORM} and/or ${BOLD}-z${NORM} are used and script is running in interactive mode, those values will be used as defaults, but may be later overwritten."
      echo
      echo -e "\tOptions required for running in non-interactive mode:"
      echo -e "\t-c\tPlastom reference sequence in FASTA format"
      echo -e "\t-x\tInput file in TSV format (output of Geneious assembly)"
      echo -e "\t-z\tInput file in FASTA format (output of Geneious assembly)"
      echo
      echo -e "\tOther optional arguments (if not provided, default values are used):"
      echo -e "\t-b\tBait length"
      echo -e "\t\tDefault value: 120 (optimal length for phylogeny, use integer between 120 and 200)"
      echo -e "\t-d\tSequence similarity between the developed probe sequences (parameter \"-c\" of cd-hit-est, see its manual for details)"
      echo -e "\t\tDefault value: 0.9 (use decimal number ranging from 0.85 to 0.95)"
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
    c)
      REFERENCECP=$OPTARG
      echo "Plastom reference: $REFERENCECP"
      ;;
    x)
      TSVLIST=$OPTARG
      echo "Input file: $TSVLIST"
      ;;
    z)
      SEQUENCES=$OPTARG
      echo "Input file: $SEQUENCES"
      ;;
    b)
      BAITL=$OPTARG
      if [[ "$BAITL" =~ ^[0-9]+$ ]] && [ "$BAITL" -ge 120 -a "$BAITL" -le 200 ]; then
	echo "Bait length: $BAITL"
	else
	  echo
	  echo "Error! For parameter \"-b\" you did not provide an integer of range from 120 to 200!"
	  echo
	  exit 1
	fi
      ;;
    d)
      CDHITSIM=$OPTARG
      if [ "$(echo 0.85 '<=' $CDHITSIM | bc -l)" = 1 ] && [ "$(echo $CDHITSIM '<=' 0.95 | bc -l)" = 1 ]; then
	echo "Sequence similarity: $CDHITSIM"
      else
	echo
	echo "Error! For parameter \"-d\" you did not provide decimal number ranging from 0.85 to 0.95!"
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

# Check if awk is available
checktools awk

# Check if sed is available
checktools sed

# Check if join is available
checktools join

# Check if cat is available
checktools cat

# Check if BLAT is available
checkblat

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

# Variables

# Function to check and read input files
# Parameters: 1) parameter for particular file; 2) name (description) of input file; 3) variable for particular file (written into $CHECKFILEREADOUT)

# Plastom reference in FASTA
CHECKFILEREADOUT=""
readinputfile -c "" $TSVLIST
TSVLIST=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

readinputfile -x "Geneious output file - TSV" $TSVLIST
TSVLIST=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Geneious output files are infiles here - consensus and unused sequences
readinputfile -z "Geneious output file - FASTA" $SEQUENCES
SEQUENCES=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Sequences converted from FASTA to tabular format - temporary file - will be deleted
SEQUENCESTAB="${SEQUENCES%.*}.tab"
# Assembled sequences in TSV - temporary file - will be deleted
SEQUENCESTABASSE="${SEQUENCES%.*}_assembled.tab"
# Unassembled sequences in TSV - temporary file - will be deleted
SEQUENCESTABUNAS="${SEQUENCESTAB%.*}_unassembled.tab"
# Filtered probes - temporary file - will be deleted
SEQUENCESPROBES600="${SEQUENCES%.*}_probes_120-600bp.tab"
# Numbers of usable contigs for joining - temporary file - will be deleted
SEQUENCESPROBES600FORJOIN="${SEQUENCESPROBES600%.*}_probes_120-600bp_fin_for_join"
# All exons ≥120 bp - temporary file - will be deleted
SEQUENCESTABASSE120="${SEQUENCES%.*}_120bp_assembled_less_than_1kb_transcript_fin.tab"
# Sorted exons ≥120 bp - temporary file - will be deleted
SEQUENCESTABASSE120SORT="${SEQUENCES%.*}_120bp_assembled_less_than_1kb_transcript_fin_sorted.tab"
# Exons ≥120 bp and all assemblies making up genes of ≥600 bp - temporary file - will be deleted
SEQUENCESPROBES120600FIN="${SEQUENCES%.*}_probes_120-600bp_fin.tab"
# Temporal files when converting from TAB to FASTA - temporary file - will be deleted
SEQUENCESPROBES120600MODIF="${SEQUENCES%.*}_probes_120-600bp_modified_fin.tab"
SEQUENCESPROBES120600ASSEM="${SEQUENCES%.*}_probes_120-600bp_assembled_fin.fasta"
SEQUENCESPROBES120600CONTIG="${SEQUENCES%.*}_probes_120-600bp_contig_fin.fasta"
# Preliminary probe sequences
PROBEPRELIM="${SEQUENCES%.*}_prelim_probe_seq.fasta"
# Sequence similarity checked by CD-HIT - temporary file - will be deleted
PROBEPRELIMCDHIT="${SEQUENCES%.*}_similarity_test"
# Assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp
PROBEPRELIMCDHIT2="${SEQUENCES%.*}_similarity_test2"
# Extracted assemblies making up genes of ≥600 bp - temporary file - will be deleted
PROBEPRELIMFORJOIN="${SEQUENCES%.*}_similarity_test_assemblies_for_join"
# Modified and sorted assemblies - temporary file - will be deleted
PROBEPRELIMSORT="${SEQUENCES%.*}_similarity_test_assemblies_sort.tab"
# All exons ≥120 bp and all assemblies making up genes of ≥600 bp - temporary file - will be deleted
PROBEPRELIMFIN="${SEQUENCES%.*}_similarity_test_assemblies_fin.tab"
# Probes in FASTA
PROBESEQUENCES="${SEQUENCES%.*}_target_enrichment_probe_sequences.fasta"

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
    $SCRIPTDIR/geneious_column_separator.pl $TSVLIST || {
      echo
      echo "${BOLD}Error!${NORM} Extraction failed. Aborting."
      echo "Please, do it manually. Required columns in $TSVLIST are \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\"."
      echo
      exit 1
      }
    TSVLIST2="${TSVLIST%.*}.columns.tsv"
    echo "File with extracted columns was saved as $TSVLIST2 for possible later usage"
  fi
echo
# { echo && echo "${BOLD}Error!${NORM} XXX failed. Aborting." && echo && exit 1; }
# Check the statistics
# Check total number of bp
echo "Total number of base pairs:"
{ cut -f6 $TSVLIST2 | awk '$1>119' | awk '{s+=$1}END{print s}'; } || { echo && echo "${BOLD}Error!${NORM} Checking statistics failed. Aborting." && echo && exit 1; }
# Check number of contigs
echo "Number of contigs:"
{ cut -f6 $TSVLIST2 | awk '$1>119' | wc -l; } || { echo && echo "${BOLD}Error!${NORM} Checking number of contigs failed. Aborting." && echo && exit 1; }
echo

# Convert FASTA to TSV
echo "Converting FASTA to TAB"
fasta2tab $SEQUENCES $SEQUENCESTAB || { echo && echo "${BOLD}Error!${NORM} Conversion of FASTA into TAB failed. Aborting." && echo && exit 1; }
echo

# Separate the assembled sequences
echo "Separating assembled sequences"
grep 'Contig' $SEQUENCESTAB > $SEQUENCESTABASSE
echo
# Separate the unassembled sequences
echo "Separating unassembled sequences"
grep -v 'Contig' $SEQUENCESTAB > $SEQUENCESTABUNAS
echo

# Filter the file with the assembled sequences – count the assemblies (the ones indicated with "Contig") making up genes of ≥960 bp / ≥600 bp, comprised of putative exons ≥120 bp
echo "Counting assembled sequences:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>959' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>959' | wc -l
echo
# Genes of ≥960 bp (exons ≥120 bp), total bp
echo "Genes of ≥960 bp (exons ≥120 bp), total bp:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
echo
# Genes of ≥600 bp (exons ≥120 bp), total bp
echo "Genes of ≥600 bp (exons ≥120 bp), total bp:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed s/_/\\t/g | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' > $SEQUENCESPROBES600
echo
# Filter the file with the unassembled sequences
echo "Filtering the file with the unassembled sequences:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABUNAS | sed s/_/\\t/g | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>399' | wc -l
echo
# Unassembled sequences making up genes of ≥400 bp
echo "Unassembled sequences making up genes of ≥400 bp:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABUNAS | sed s/_/\\t/g | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
awk '{print $1"\t"length($2)}' $SEQUENCESTABUNAS | sed s/_/\\t/g | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
echo

# Part 4: Create the final .fasta file for the Hyb-Seq probes.

# Extract and sort the assemblies making up genes of ≥600 bp
echo "Extracting and sorting the assemblies making up genes of ≥600 bp"
sed s/^/Assembly_/ $SEQUENCESPROBES600 | cut -f1 -d' ' | sort -k1,1 > $SEQUENCESPROBES600FORJOIN
echo
# Make a file with all exons ≥120 bp
echo "Selecting ≥120 bp exons"
awk '{print $1"\t"length($2)"\t"$2}' $SEQUENCESTABASSE | awk '$2>119' > $SEQUENCESTABASSE120
echo
# Make the assembly number the first field and sort
echo "Sorting exons ≥120 bp"
sed 's/^.*\(Assembly\)/\1/' $SEQUENCESTABASSE120 | sed s/_C/\\tC/ | sort -k1,1 > $SEQUENCESTABASSE120SORT
echo
# Make a file with all exons ≥120 bp and all assemblies making up genes of ≥600 bp
echo "Selecting all exons ≥120 bp and all assemblies making up genes of ≥600 bp"
join $SEQUENCESPROBES600FORJOIN $SEQUENCESTABASSE120SORT > $SEQUENCESPROBES120600FIN
echo
# Convert .tab to .fasta
echo "Converting TAB to FASTA"
sed 's/ /_/' $SEQUENCESPROBES120600FIN | sed 's/ /_/' > $SEQUENCESPROBES120600MODIF
sed 's/^/>/' $SEQUENCESPROBES120600MODIF | tr " " "\n" > $SEQUENCESPROBES120600ASSEM
# Remaining assemblies have to be selected and added to the .fasta file of the probes:
grep -v Contig $SEQUENCESTABASSE120 | awk '$2>599' | sed 's/^/>/' | sed 's/\\t/_/' | tr " " "\n" > $SEQUENCESPROBES120600CONTIG
echo
# Combine the two .fasta files:
echo "Writing FASTA file with preliminary probe sequences"
cat $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG > $PROBEPRELIM
echo
echo "Preliminary probe sequences saved as $PROBEPRELIM for possible later usage"

# Part 5: Make the final quality control of the probe sequences before sending them to company for bait synthesis

# Check for sequence similarity between the developed probe sequences with CD-HIT-EST
echo "Checking sequence similarity between the developed probe sequences"
cd-hit-est -i $PROBEPRELIM -o $PROBEPRELIMCDHIT -c $CDHITSIM
echo
# One of the three outfiles is a FASTA file, it has to be converted to TAB
echo "Converting FASTA to TAB"
fasta2tab $PROBEPRELIMCDHIT $PROBEPRELIMCDHIT.txt || { echo && echo "${BOLD}Error!${NORM} Conversion of FASTA into TAB failed. Aborting." && echo && exit 1; }
echo
# Count all assemblies, comprised of putative exons ≥120 bp
echo "Counting all assemblies, comprised of putative exons ≥120 bp:"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | awk '{s+=$2;c++}END{print s}'
echo
# Count the assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp
echo "Counting the assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp:"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | sed s/_/\\t/g | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | sed s/_/\\t/g | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
echo
echo "Writing those assemblies into temporal file"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | sed s/_/\\t/g | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' > $PROBEPRELIMCDHIT2
echo
# Extract and sort the assemblies making up genes of ≥600 bp
echo "Extract and sort the assemblies making up genes of ≥600 bp:"
sed s/^/Assembly_/ $PROBEPRELIMCDHIT2 | cut -f1 -d' ' | sort -k1,1 > $PROBEPRELIMFORJOIN
echo
# Modify the assembly number and sort
echo "Modify the assembly number and sort"
sed s/_C/\\tC/ $PROBEPRELIMCDHIT.txt | sort -k1,1 > $PROBEPRELIMSORT
echo
# Make a file with all exons ≥120 bp and all assemblies making up genes of ≥600 bp
echo "Joining all exons ≥120 bp and all assemblies making up genes of ≥600 bp"
join $PROBEPRELIMFORJOIN $PROBEPRELIMSORT > $PROBEPRELIMFIN
echo
# Convert TAB to FASTA
echo "Converting TAB to FASTA"
sed s/' C'/_/ $PROBEPRELIMFIN | sed s/ontig/Contig/ | sed s/^/'>'/ | sed s/' '/\\n/ > $PROBESEQUENCES
echo
# Calculating of the total number of base pairs
echo "Calculating of the total number of base pairs"
echo "Converting FASTA to TAB"
fasta2tab $PROBESEQUENCES $PROBESEQUENCESNUM || { echo && echo "${BOLD}Error!${NORM} Conversion of FASTA into TAB failed. Aborting." && echo && exit 1; }
echo
echo "Total number of base pairs:"
awk '{print $1"\t"length($2)}' $PROBESEQUENCESNUM | awk '{s+=$2;c++}END{print s}'
echo
# Remove remaining cp/mt genes from probe set
echo "Removing remaining cp/mt genes from probe set"
blat -t=dna -q=dna -out=pslx $REFERENCECP $PROBESEQUENCES $PROBESEQUENCES.target_enrichment_probe_sequences_final.pslx
echo

# Remove temporal files
echo "Removing unneeded temporal files"
rm $SEQUENCESTAB $SEQUENCESTABASSE $SEQUENCESTABUNAS $SEQUENCESPROBES600 $SEQUENCESPROBES600FORJOIN $SEQUENCESTABASSE120 $SEQUENCESTABASSE120SORT $SEQUENCESPROBES120600FIN $SEQUENCESPROBES120600MODIF $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG $PROBEPRELIMCDHIT $PROBEPRELIMFORJOIN $PROBEPRELIMSORT $PROBEPRELIMFIN $PROBESEQUENCESNUM
echo

echo "Success!"
echo
echo "Final output file was written as $PROBESEQUENCES.target_enrichment_probe_sequences_final.pslx"
echo
echo "This file contains the probe sequences"
echo

echo "Script exited successfully..."
echo

exit
