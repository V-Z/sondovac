#!/bin/bash

# Determine script's directory
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load functions shared by both parts, introductory message
source $SCRIPTDIR/sondovac_functions || {
  echo
  echo "Fatal error!"
  echo "Unable to load file \"sondovac_functions\" with required functions!"
  echo "It must be in same directory as \"$0\""
  echo "Check it and if needed download again whole script from"
  echo "https://github.com/V-Z/sondovac/"
  echo
  exit 1
  }

echo "This is part B of the pipeline."
echo
echo "This part processes assembly output from Geneious and produces the final list of"
echo "low-copy nuclear probe sequences."

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
      generaloptions
      echo
      echo -e "\tIf options ${BOLD}-c${NORM}, ${BOLD}-x${NORM} and/or ${BOLD}-z${NORM} are used and script is running in"
      echo -e "\t  interactive mode, those values will be used as defaults, but may be"
      echo -e "\t  later overwritten."
      echo
      echo -e "\tOptions required for running in non-interactive mode:"
      echo -e "\t-c\tPlastome reference sequence input file in FASTA format."
      echo -e "\t-x\tInput file in TSV format (output of Geneious assembly)."
      echo -e "\t-z\tInput file in FASTA format (output of Geneious assembly)."
      echo
      echo -e "\tOther optional arguments (if not provided, default values are used):"
      echo -e "\t-b\tBait length"
      echo -e "\t\tDefault value: 120 (optimal length for phylogeny, use integer"
      echo -e "\t\t  between 120 and 200)."
      echo -e "\t-d\tSequence similarity between the developed probe sequences"
      echo -e "\t\t  (parameter \"-c\" of cd-hit-est, see its manual for details)."
      echo -e "\t\tDefault value: 0.9 (use decimal number ranging from 0.85 to 0.95)."
      echo -e "\t${BOLD}WARNING!${NORM} If parameters ${BOLD}-b${NORM} or ${BOLD}-d${NORM} are not provided, default values are"
      echo -e "\t  taken and it is not possible to change them later (not even in"
      echo -e "\t  interactive mode)."
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
	  echo "${BOLD}Error!${NORM} For parameter \"-b\" you did not provide an integer ranging from 120 to 200!"
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
	echo "${BOLD}Error!${NORM} For parameter \"-d\" you did not provide decimal number ranging from 0.85"
	echo "  to 0.95!"
	echo
	exit 1
      fi
      ;;
    ?)
      echo
      echo "Invalid option(s)!"
      echo "See \"$0 -h\" for usage options."
      echo
      exit 1
      ;;
    esac
  done

# Check if user didn't use together -n and -i
checkmodef

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

# Check if cut is available
checktools cut

# Check if awk is available
checktools awk

# Check if wc is available
checktools wc

# Check if sed is available
checktools sed

# Check if sort is available
checktools sort

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
  echo "Compiling CD-HIT from source code..." &&
  echo &&
  cd $1 &&
  make -s openmp=yes || { echo "There is no MPI available - no multi-thread support." && make -s; } &&
  cp -a *.pl $BIN/ &&
  cp -a cd-hit* $BIN/ &&
  cd $WORKDIR &&
  echo &&
  echo "\"CD-HIT\" is available. OK."
  } || {
    echo
    echo "Compilation failed. Please go to https://github.com/weizhongli/cdhit"
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
  echo "Type \"C\" to compile CD-HIT 4.6.4 from source available together with this script."
  echo "Type \"S\" to download latest CD-HIT source from"
  echo "  https://github.com/weizhongli/cdhit and compile it"
  echo "Type \"B\" to copy CD-HIT 4.6.4 binary available together with the script"
  echo "  (recommended, available for Linux and Mac OS X)."
  echo "Type \"M\" for manual installation - script will exit and you will have to install"
  echo "  CD-HIT yourselves."
  read CDHIT
  while :
  do
    case "$CDHIT" in
      C|c)
	compilecdhit $SCRIPTDIR/src/cd-hit-v4.6.4-2015-0603
	break
	;;
      S|s)
	downloaderselector
	checktools unzip
	$DOWNLOADER cd-hit-master.zip https://github.com/weizhongli/cdhit/archive/master.zip
	unzip -nq cd-hit-master.zip
	compilecdhit cdhit-master
	break
	;;
      B|b)
	echo "Copying CD-HIT binaries"
	case "$OS" in
	  Mac)
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/cd-hit* $BIN/
	    cp -p $SCRIPTDIR/pkgs/macosx/bin/*.pl $BIN/
	    ;;
	  Linux)
	    cp -p $SCRIPTDIR/pkgs/linux64b/bin/cd-hit* $BIN/
	    cp -p $SCRIPTDIR/pkgs/linux64b/bin/*.pl $BIN/
	    ;;
	  *) echo
	    echo "Binary is not available for $OS $OSB."
	    echo
	    compilecdhit $SCRIPTDIR/src/cd-hit-v4.6.4-2015-0603
	    ;;
	esac
	break
	;;
      M|m)
	echo
	echo "Please, go to http://weizhongli-lab.org/cd-hit/ and install CD-HIT and ensure"
	echo "  it is in PATH."
	echo
	exit
	;;
      *) echo "Wrong option. Use C, S, B or M." && read CDHIT;;
    esac
  done
  }

# Input files
CHECKFILEREADOUT=""

# Plastom reference in FASTA
readinputfile -c "plastome reference sequence input file in FASTA format" $TSVLIST
REFERENCECP=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Geneious output files are infiles here - consensus and unused sequences (TSV)
readinputfile -x "input file in TSV format (output of Geneious assembly)" $TSVLIST
TSVLIST=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Geneious output files are infiles here - consensus and unused sequences (FASTA)
readinputfile -z "input file in FASTA format (output of Geneious assembly)" $SEQUENCES
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
# Preliminary probe sequences - temporary file - will be deleted
PROBEPRELIM0="${SEQUENCES%.*}_prelim_probe_seq0.fasta"
# Preliminary probe sequences - corrected labels
PROBEPRELIM="${SEQUENCES%.*}_prelim_probe_seq.fasta"
# Sequence similarity checked by CD-HIT - temporary file - will be deleted
PROBEPRELIMCDHIT="${SEQUENCES%.*}_sim_test.fasta"
# Assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp
PROBEPRELIMCDHIT2="${SEQUENCES%.*}_similarity_test.fasta"
# Extracted assemblies making up genes of ≥600 bp - temporary file - will be deleted
PROBEPRELIMFORJOIN="${SEQUENCES%.*}_similarity_test_assemblies_for_join"
# Modified and sorted assemblies - temporary file - will be deleted
PROBEPRELIMSORT="${SEQUENCES%.*}_similarity_test_assemblies_sort.tab"
# All exons ≥120 bp and all assemblies making up genes of ≥600 bp - temporary file - will be deleted
PROBEPRELIMFIN="${SEQUENCES%.*}_similarity_test_assemblies_fin.tab"
# Probes in FASTA
PROBESEQUENCES="${SEQUENCES%.*}_target_enrichment_probe_sequences.fasta"
# Probes in tab - temporary file - will be deleted
PROBESEQUENCESNUM="${SEQUENCES%.*}_target_enrichment_probe_sequences.tab"
# Possible cpDNA genes in probe set
PROBESEQUENCESCP="${SEQUENCES%.*}_possible_cp_dna_genes_in_probe_set.pslx"

# Assemble the obtained sequences in contigs

# Step 8: Retention of those contigs that comprise exons ≥ bait length (default is 120 bp) and have a certain locus length

echo "Step 8 of the pipeline - retention of those contigs that comprise exons ≥ bait"
echo "length ($BAITL bp) and have a certain locus length"
echo

# Check if TSV output of Geneious contains at least requested columns
REQUIREDCOLS="Required columns in `echo $TSVLIST` are \"# Sequences\",\n  \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\"\n  and \"Sequence Length\". Please, export the TSV file again."
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
	      echo "${BOLD}Error!${NORM} Column \"Sequence Length\" is missing!"
	      echo -e "$REQUIREDCOLS"
	      exit 1
	    fi
	  else
	    echo "${BOLD}Error!${NORM} Column \"Name\" is missing!"
	    echo -e "$REQUIREDCOLS"
	    exit 1
	  fi
	else
	  echo "${BOLD}Error!${NORM} Column \"Mean Coverage\" is missing!"
	  echo -e "$REQUIREDCOLS"
	  exit 1
	fi
      else
	echo "${BOLD}Error!${NORM} Column \"Description\" is missing!"
	echo -e "$REQUIREDCOLS"
	exit 1
      fi
    else
      echo "${BOLD}Error!${NORM} Column \"Pairwise Identity\" is missing!"
      echo -e "$REQUIREDCOLS"
      exit 1
    fi
  else
    echo "${BOLD}Error!${NORM} Column \"# Sequences\" is missing!"
    echo -e "$REQUIREDCOLS"
    exit 1
  fi
echo

# Check if TSV output of Geneious contains only required columns or more
if egrep -q "# Sequences[[:blank:]]+% Pairwise Identity[[:blank:]]+Description[[:blank:]]+Mean Coverage[[:blank:]]+Name[[:blank:]]+Sequence Length" $TSVLIST
  then
    echo "OK, $TSVLIST is correct input file."
    TSVLIST2=$TSVLIST
  else
    echo "Input file $TSVLIST seems to contain more columns than required."
    echo "Needed columns will be extracted."
    $SCRIPTDIR/geneious_column_separator.pl $TSVLIST || {
      echo
      echo "${BOLD}Error!${NORM} Extraction failed. Aborting."
      echo "Either script $SCRIPTDIR/geneious_column_separator.pl"
      echo "  is missing or there is something wrong with $TSVLIST"
      echo "Please, prepare required file manually."
      echo -e "$REQUIREDCOLS"
      echo
      exit 1
      }
    TSVLIST2="${TSVLIST%.*}.columns.tsv"
    echo "File with extracted columns was saved as"
    echo "$TSVLIST2 for possible later usage."
    confirmgo
  fi
echo

# Check the statistics

echo "Assembly statistics"
echo

# Check total number of bp
echo "Total number of base pairs:"
{ cut -f6 $TSVLIST2 | awk '$1>119' | awk '{s+=$1}END{print s}'; } || {
  echo
  echo "${BOLD}Error!${NORM} Checking statistics failed. Aborting. Check if file"
  echo "$TSVLIST2 is correct TSV file containing all required columns:"
  echo -e "$REQUIREDCOLS"
  echo
  exit 1
  }
confirmgo
echo

# Check number of contigs
echo "Number of contigs:"
{ cut -f6 $TSVLIST2 | awk '$1>119' | wc -l; } || {
  echo
  echo "${BOLD}Error!${NORM} Checking number of contigs failed. Aborting. Check if file"
  echo "$TSVLIST2 is correct TSV file containing all required columns"
  echo -e "$REQUIREDCOLS"
  echo
  exit 1
  }
confirmgo
echo

# Convert FASTA to TSV
echo "Converting FASTA to TAB"
fasta2tab $SEQUENCES $SEQUENCESTAB || {
  echo
  echo "${BOLD}Error!${NORM} Conversion of FASTA into TAB failed. Aborting."
  echo "Check if file $SEQUENCES is correct FASTA file."
  echo
  exit 1
  }
echo

# Separate the assembled sequences
echo "Separating assembled sequences"
grep '[Aa]ssembly' $SEQUENCESTAB > $SEQUENCESTABASSE
echo

# Separate the unassembled sequences
echo "Separating unassembled sequences"
grep -v '[Aa]ssembly' $SEQUENCESTAB > $SEQUENCESTABUNAS
echo

# Filter the file with the assembled sequences – count the assemblies (the ones indicated with "Contig") making up genes of ≥960 bp / ≥600 bp, comprised of putative exons ≥120 bp
echo "Number of assembled sequences:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/_/\t/g' | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>959' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/_/\t/g' | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>959' | wc -l
confirmgo
echo

# Genes of ≥960 bp (exons ≥120 bp), total bp
echo "Genes of ≥960 bp (exons ≥120 bp), total bp:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/_/\t/g' | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/_/\t/g' | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
confirmgo
echo

# Genes of ≥600 bp (exons ≥120 bp), total bp
echo "Genes of ≥600 bp (exons ≥120 bp), total bp:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/_/\t/g' | cut -f6,9 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' > $SEQUENCESPROBES600
confirmgo
echo

# Filter the file with the unassembled sequences
echo "Filtering the file with the unassembled sequences:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABUNAS | sed 's/_/\t/g' | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>399' | wc -l
confirmgo
echo

# Unassembled sequences making up genes of ≥400 bp
echo "Unassembled sequences making up genes of ≥400 bp:"
awk '{print $1"\t"length($2)}' $SEQUENCESTABUNAS | sed 's/_/\t/g' | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
awk '{print $1"\t"length($2)}' $SEQUENCESTABUNAS | sed 's/_/\t/g' | cut -f1,4 | awk '$2>119' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
confirmgo
echo

# Create the final FASTA file for the Hyb-Seq probes

# Extract and sort the assemblies making up genes of ≥600 bp
echo "Extracting and sorting the assemblies making up genes of ≥600 bp"
sed 's/^/Assembly_/' $SEQUENCESPROBES600 | cut -f1 -d " " | sort -k1,1 > $SEQUENCESPROBES600FORJOIN
echo

# Make a file with all exons ≥120 bp
echo "Selecting ≥120 bp exons"
awk '{print $1"\t"length($2)"\t"$2}' $SEQUENCESTABASSE | awk '$2>119' > $SEQUENCESTABASSE120
echo

# Make the assembly number the first field and sort
echo "Sorting exons ≥120 bp"
sed 's/^.*\(Assembly\)/\1/' $SEQUENCESTABASSE120 | sed 's/_C/\tC/' | sort -k1,1 > $SEQUENCESTABASSE120SORT
echo

# Make a file with all exons ≥120 bp and all assemblies making up genes of ≥600 bp
echo "Selecting all exons ≥120 bp and all assemblies making up genes of ≥600 bp"
join $SEQUENCESPROBES600FORJOIN $SEQUENCESTABASSE120SORT > $SEQUENCESPROBES120600FIN
echo

# Convert TAB to FASTA
echo "Converting TAB to FASTA"
sed 's/ /_/' $SEQUENCESPROBES120600FIN | sed 's/ /_/' > $SEQUENCESPROBES120600MODIF
sed 's/^/>/' $SEQUENCESPROBES120600MODIF | sed 's/ /\n/' > $SEQUENCESPROBES120600ASSEM

# Remaining assemblies have to be selected and added to the .fasta file of the probes:
grep -v '[Cc]ontig' $SEQUENCESTABASSE120 | awk '$2>599' | sed 's/^/>/' | sed 's/\t/_/' | sed 's/\t/\n/' > $SEQUENCESPROBES120600CONTIG
echo

# Combine the two FASTA files
echo "Writing FASTA file with preliminary probe sequences"
cat $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG > $PROBEPRELIM0

# Ensure all sequences have correct labeling
echo "Ensuring all sequences have correct labels"
sed 's/^>.\+\(Assembly_[0-9]\+_\)/>\1Contig_0_/' $PROBEPRELIM0 > $PROBEPRELIM
echo
echo "Preliminary probe sequences saved as $PROBEPRELIM for possible later usage"
confirmgo

# Step 9: Make the final quality control of the probe sequences before sending them to company for bait synthesis

echo
echo "Step 9 of the pipeline - removal of probe sequences sharing ≥90% sequence"
echo "similarity"
echo

# Check for sequence similarity between the developed probe sequences with CD-HIT-EST
echo "Checking sequence similarity between the developed probe sequences"
cd-hit-est -i $PROBEPRELIM -o $PROBEPRELIMCDHIT -c $CDHITSIM

# Step 10: Retention of those contigs that comprise exons ≥ bait length (default is 120 bp) and have a certain locus length

echo
echo "Step 10 of the pipeline - retention of those contigs that comprise exons ≥ bait"
echo "length ($BAITL bp) and have a certain locus length"
echo

# One of the three outfiles is a FASTA file, it has to be converted to TAB
echo "Converting FASTA to TAB"
fasta2tab $PROBEPRELIMCDHIT $PROBEPRELIMCDHIT.txt || {
  echo
  echo "${BOLD}Error!${NORM} Conversion of FASTA into TAB failed. Aborting."
  echo "Check if file $PROBEPRELIMCDHIT is correct FASTA file."
  echo
  exit 1
  }
echo

# Count all assemblies, comprised of putative exons ≥120 bp
echo "Number of all assemblies, comprised of putative exons ≥120 bp:"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | awk '{s+=$2;c++}END{print s}'
echo
confirmgo

# Count the assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp
echo "Number of the assemblies making up genes of ≥600 bp, comprised of putative exons ≥120 bp:"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | sed 's/_/\t/g' | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | awk '{s+=$3;c++}END{print s}'
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | sed 's/_/\t/g' | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' | wc -l
echo
confirmgo

echo "Writing the assemblies into temporal file"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT.txt | sed 's/_/\t/g' | cut -f2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>599' > $PROBEPRELIMCDHIT2
echo

# Extract and sort the assemblies making up genes of ≥600 bp
echo "Extract and sort the assemblies making up genes of ≥600 bp:"
sed 's/^/Assembly_/' $PROBEPRELIMCDHIT2 | cut -f1 -d " " | sort -k1,1 > $PROBEPRELIMFORJOIN
echo

# Modify the assembly number and sort
echo "Modify the assembly number and sort"
sed 's/_C/\tC/' $PROBEPRELIMCDHIT.txt | sort -k1,1 > $PROBEPRELIMSORT
echo

# Make a file with all exons ≥120 bp and all assemblies making up genes of ≥600 bp
echo "Joining all exons ≥120 bp and all assemblies making up genes of ≥600 bp"
join $PROBEPRELIMFORJOIN $PROBEPRELIMSORT > $PROBEPRELIMFIN
echo

# Convert TAB to FASTA
echo "Converting TAB to FASTA"
sed 's/^\(.\+\) \(Contig\)/>\1_\2/' $PROBEPRELIMFIN | sed 's/ /\n/' > $PROBESEQUENCES
echo

# Calculating of the total number of base pairs
echo "Calculating of the total number of base pairs"
echo "Converting FASTA to TAB"
fasta2tab $PROBESEQUENCES $PROBESEQUENCESNUM || {
  echo
  echo "${BOLD}Error!${NORM} Conversion of FASTA into TAB failed. Aborting."
  echo
  exit 1
  }

echo
echo "Total number of base pairs:"
awk '{print $1"\t"length($2)}' $PROBESEQUENCESNUM | awk '{s+=$2;c++}END{print s}'
confirmgo
echo

echo "${BOLD}Success!${NORM}"
echo
echo "================================================================================"
echo "Final output file was written as"
echo "${BOLD}$PROBESEQUENCES${NORM}"
echo "This file contains the probe sequences."
echo "================================================================================"
confirmgo

# Step 11 - removal of possible cpDNA sequences in final probe list

echo
echo "Step 11 of the pipeline - removal of probe sequences sharing ≥90% sequence"
echo "similarity with the plastome reference"
echo

# Remove remaining cp genes from probe set
echo "Removing remaining plastid genes from probe set"
blat -t=dna -q=dna -out=pslx $REFERENCECP $PROBESEQUENCES $PROBESEQUENCESCP
echo

echo "================================================================================"
echo "File $PROBESEQUENCESCP"
echo "contains possible plastid genes in final probe set."
echo "We recommend to remove those genes from final probe set in file"
echo "$PROBESEQUENCES."
echo "================================================================================"
confirmgo
echo

# Remove temporal files
echo "Removing unneeded temporal files"
rm $SEQUENCESTAB $SEQUENCESTABASSE $SEQUENCESTABUNAS $SEQUENCESPROBES600 $SEQUENCESPROBES600FORJOIN $SEQUENCESTABASSE120 $SEQUENCESTABASSE120SORT $SEQUENCESPROBES120600FIN $SEQUENCESPROBES120600MODIF $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG $PROBEPRELIMCDHIT $PROBEPRELIMFORJOIN $PROBEPRELIMSORT $PROBEPRELIMFIN $PROBESEQUENCESNUM || {
  echo
  echo "${BOLD}Error!${NORM} Removal of temporal files failed. Remove following files manually:"
  echo "$SEQUENCESTAB, $SEQUENCESTABASSE, $SEQUENCESTABUNAS,"
  echo "$SEQUENCESPROBES600, $SEQUENCESPROBES600FORJOIN, $SEQUENCESTABASSE120,"
  echo "$SEQUENCESTABASSE120SORT, $SEQUENCESPROBES120600FIN, $SEQUENCESPROBES120600MODIF,"
  echo "$SEQUENCESPROBES120600ASSEM, $SEQUENCESPROBES120600CONTIG, $PROBEPRELIMCDHIT,"
  echo "$PROBEPRELIMFORJOIN, $PROBEPRELIMSORT,"
  echo "$PROBEPRELIMFIN and $PROBESEQUENCESNUM."
  confirmgo
  }

# List kept files which user can use for another analysis
echo
echo "Following files are kept for possible later usage (see manual for details):"
echo "================================================================================"
echo "1) Preliminary probe sequences:"
echo "$PROBEPRELIM"
echo "2) Contigs comprising exons ≥ bait length and have a certain total locus length:"
echo "$PROBEPRELIMCDHIT2"
echo "3) Possible plastid genes in final probe set:"
echo "$PROBESEQUENCESCP"
echo "4) Final probe sequences in FASTA format:"
echo "$PROBESEQUENCES"
echo "================================================================================"
confirmgo

echo
echo "Script exited successfully..."
echo

exit
