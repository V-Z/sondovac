#!/bin/bash

# Determine script's directory
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load aliases to replace Mac OS X outdated tools by those installed by Homebrew
shopt -s expand_aliases
unalias -a
source $SCRIPTDIR/mac_aliases

# Load functions shared by both parts, introductory message
source $SCRIPTDIR/sondovac_functions || {
	echo
	echo "Fatal error!"
	echo "Unable to load file \"sondovac_functions\" with required functions!"
	echo "It must be in same directory as \"$0\""
	echo "Check it and, if needed, download again whole script from"
	echo "https://github.com/V-Z/sondovac/"
	echo
	exit 1
	}

echo "${REDF}This is part B of the pipeline.${NORM}"
echo
echo "This part processes assembly output from Geneious and produces the final list of"
echo "low-copy nuclear probe sequences."

# Default values
# Counter if not both -i and -n options are used
CHECKMODE=0
# If not specifying explicitly otherwise (using -n), running in interactive mode
STARTINI="I"
# Minimum bait/exon length
BAITL=120
# CD-HIT sequence similarity
CDHITSIM=0.9
# BLAT -minIdentity between the probe sequences and the plastome reference
BLATIDENT=90
# Minimum total locus length
MINLOCUSLENGTH=600
# Default name of output files
OUTPUTFILENAME="output"
# Add also unassembled sequences longer than minimal length
REMAINING="YES"

# Create empty variables for file names
REFERENCECP=""
REFERENCECP0=""
TSVLIST=""
SEQUENCES=""
SEQUENCES0=""

# Parse initial arguments
while getopts "hvulrpeo:inc:x:z:b:d:y:k:" START; do
  case "$START" in
    h|v)
      generaloptions
      echo
      echo -e "\tIf options ${BOLD}-c${NORM}, ${BOLD}-x${NORM} and/or ${BOLD}-z${NORM} are used and script is running in"
      echo -e "\t  interactive mode, those values will be used as defaults, but may be"
      echo -e "\t  later overwritten."
      echo
      echo -e "\tOptions required for running in non-interactive mode:"
      echo -e "\t${REDF}-c${NORM}\t${CYAF}Plastome reference sequence${NORM} input file in FASTA format."
      echo -e "\t${REDF}-x${NORM}\t${CYAF}Input file in TSV${NORM} format (output of Geneious assembly)."
      echo -e "\t${REDF}-z${NORM}\t${CYAF}Input file in FASTA${NORM} format (output of Geneious assembly)."
      echo
      echo -e "\tOther optional arguments (if not provided, default values are used):"
      echo -e "\t${REDF}-b${NORM}\t${CYAF}Bait length${NORM}"
      echo -e "\t\tDefault value: 120 (preferred length for phylogeny, use any of"
      echo -e "\t\t  values 80, 100 or 120)."
      echo -e "\t${REDF}-d${NORM}\t${CYAF}Sequence similarity between the developed probe sequences${NORM}"
      echo -e "\t\t  (parameter \"-c\" of cd-hit-est, see its manual for details)."
      echo -e "\t\tDefault value: 0.9 (use decimal number ranging from 0.85 to"
      echo -e "\t\t  0.95)."
      echo -e "\t${REDF}-y${NORM}\t${CYAF}Sequence similarity between the probes and plastome reference${NORM}"
      echo -e "\t\t  searching for possible plastid genes in probe set (parameter"
      echo -e "\t\t  \"-minIdentity\" of BLAT, see its manual for details)."
      echo -e "\t\tDefault value: 90 (integer ranging from 85 to 95)."
      echo -e "\t${REDF}-k${NORM}\t${CYAF}Minimum total locus length.${NORM}"
      echo -e "\t\tDefault value: 600. Allowed values are 360, 480, 600, 720, 840,"
      echo -e "\t\t  960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920 and 2040."
      echo -e "\t\t  When running in interactive mode, the user will be asked"
      echo -e "\t\t  which value to use. A table summarizing the total number"
      echo -e "\t\t  of LCN loci and the total number of base pairs for these"
      echo -e "\t\t  values will be displayed to facilitate this choice."
      echo -e "\t${BOLD}WARNING!${NORM} If parameters ${BOLD}-b${NORM}, ${BOLD}-d${NORM} or ${BOLD}-y${NORM} are not provided, default values"
      echo -e "\t  are taken, and it is not possible to change them later (not even in"
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
    o)
      OUTPUTFILENAME=$OPTARG
      echo "Output files will start name with ${REDF}$OUTPUTFILENAME${NORM}"
      ;;
    i)
      echo "${CYAF}Running in interactive mode...${NORM}"
      STARTINI="I"
      CHECKMODE=$((CHECKMODE+1))
      ;;
    n)
      echo "${CYAF}Running in non-interactive mode...${NORM}"
      STARTINI="N"
      CHECKMODE=$((CHECKMODE+1))
      ;;
    c)
      REFERENCECP0=$OPTARG
      echo "Plastome reference: ${REDF}$REFERENCECP${NORM}"
      ;;
    x)
      TSVLIST=$OPTARG
      echo "Input file: ${REDF}$TSVLIST${NORM}"
      ;;
    z)
      SEQUENCES0=$OPTARG
      echo "Input file: ${REDF}$SEQUENCES${NORM}"
      ;;
    b)
      BAITL=$OPTARG
      # Check if provided value makes sense
      case "$BAITL" in
        80) BAITL=80;;
        100) BAITL=100;;
        120) BAITL=120;;
        *) echo
        echo "${REDF}${BOLD}Error!${NORM} For parameter \"-b\" you did not provide any of values 80, 100 or 120!"
        echo
        exit 1
        esac
      echo "Bait length: ${REDF}$BAITL${NORM}"
      ;;
    d)
      CDHITSIM=$OPTARG
      # Check if provided value makes sense
      if [ "$(echo 0.85 '<=' $CDHITSIM | bc -l)" = 1 ] && [ "$(echo $CDHITSIM '<=' 0.95 | bc -l)" = 1 ]; then
        echo "Sequence similarity: $CDHITSIM"
        else
        echo
        echo "${REDF}${BOLD}Error!${NORM} For parameter \"-d\" you did not provide decimal number ranging from 0.85"
        echo "  to 0.95!"
        echo
        exit 1
        fi
      ;;
    y)
      BLATIDENT=$OPTARG
      # Check if provided value makes sense
      if [[ "$BLATIDENT" =~ ^[0-9]+$ ]] && [ "$BLATIDENT" -ge 85 -a "$BLATIDENT" -le 95 ]; then
        echo "BLAT score for similarity between unique transcripts and genome skim data: $BLATIDENT"
        else
        echo
        echo "${REDF}${BOLD}Error!${NORM} For parameter \"-y\" you did not provide an integer of range from 85 to 95!"
        echo
        exit 1
        fi
      ;;
    k)
      MINLOCUSLENGTH=$OPTARG
      # Check if provided value makes sense
      case "$MINLOCUSLENGTH" in
        360) MINLOCUSLENGTH=360;;
        480) MINLOCUSLENGTH=480;;
        600) MINLOCUSLENGTH=600;;
        720) MINLOCUSLENGTH=720;;
        840) MINLOCUSLENGTH=840;;
        960) MINLOCUSLENGTH=960;;
        1080) MINLOCUSLENGTH=1080;;
        1200) MINLOCUSLENGTH=1200;;
        1320) MINLOCUSLENGTH=1320;;
        1440) MINLOCUSLENGTH=1440;;
        1560) MINLOCUSLENGTH=1560;;
        1680) MINLOCUSLENGTH=1680;;
        1800) MINLOCUSLENGTH=1800;;
        1920) MINLOCUSLENGTH=1920;;
        2040) MINLOCUSLENGTH=2040;;
        *) echo
        echo "${REDF}${BOLD}Error!${NORM} For parameter \"-k\" you did not provide any of values 360, 480, 600, 720,"
        echo "  840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920 or 2040!"
        echo
        exit 1
        esac
      echo "Minimum exon length: ${REDF}$MINLOCUSLENGTH${NORM}"
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

# Set bait length
BAITLN=$(expr $BAITL - 1)

# Check if user didn't use together -n and -i
checkmodef

# Check which operating system the script is running on

# Ensure user reads introductory information
confirmgo

# Check operating system
oscheck

# Set variables for working directory and ${BOLD}PATH${NORM}
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

# Check if Python is available
checktools python

# Function to compile CD-HIT
function compilecdhit {
	{
	checktools make &&
	checktools g++ &&
	echo &&
	echo "Compiling \"${REDF}CD-HIT${NORM}\" from source code..." &&
	echo &&
	cd $1 &&
	make -s openmp=yes || { echo "${CYAF}There is no MPI available${NORM} - no multi-thread support." && make -s openmp=no; } &&
	cp -a *.pl $BIN/ || echo "No Perl scripts in this build..." &&
	cp -a cd-hit* $BIN/ &&
	cd $WORKDIR &&
	echo &&
	echo "\"CD-HIT\" is available. ${GREF}OK.${NORM}"
		} || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please go to ${REDF}https://github.com/weizhongli/cdhit/releases${NORM}"
		echo "download cd-hit-*.tgz, compile it and ensure it is in ${BOLD}PATH${NORM}."
		echo "Check last error messages to find out why compilation failed."
		echo
		exit 1
		}
	}

# Check if cd-hit-est is available
{ command -v cd-hit-est >/dev/null 2>&1 && echo "\"${REDF}cd-hit-est${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
	echo "${CYAF}CD-HIT is required but not installed or available in ${BOLD}PATH${NORM}.${NORM}"
	if [ "$STARTINI" == "I" ]; then
		echo
		echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile \"CD-HIT\" 4.6.4 from source available together with this script.${NORM}"
		echo "Type \"${REDF}S${NORM}\" ${CYAF}to download latest \"CD-HIT\" source${NORM} from"
		echo "  ${REDF}https://github.com/weizhongli/cdhit${NORM} and compile it"
		echo "Type \"${REDF}B${NORM}${NORM}\" ${CYAF}to copy \"CD-HIT\" 4.6.4 binary${NORM} available together with the script"
		echo "  (recommended, available for Linux and Mac OS X)."
		echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
		echo "  See \"${REDF}brew info homebrew/science/cd-hit${NORM}\" for more details."
		echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit, and you will have to install"
		echo "  \"CD-HIT\" yourself."
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
						$DOWNLOADER cd-hit-master.zip https://github.com/weizhongli/cdhit/archive/master.zip || {
							echo
							echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to"
							echo "  ${REDF}https://github.com/weizhongli/cdhit/releases/${NORM},"
							echo "  download latest cd-hit-*.tar.gz and compile it manually."
							echo
							exit 1
							}
						unzip -nq cd-hit-master.zip
						compilecdhit cdhit-master
						break
						;;
					B|b)
						echo "Copying \"${REDF}CD-HIT${NORM}\" binaries"
						case "$OS" in
							Mac)
								cp -pr $SCRIPTDIR/pkgs/macosx/bin/cd-hit* $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/*.pl $BIN/ || echo
								;;
							Linux)
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/cd-hit* $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/*.pl $BIN/
								;;
							*) echo
								echo "Binary is not available for ${REDF}$OS $OSB${NORM}."
								echo
								compilecdhit $SCRIPTDIR/src/cd-hit-v4.6.4-2015-0603
								;;
							esac
						break
						;;
					H|h)
						if [ "$OS" == "Mac" ]; then
							{ echo "Installing \"${REDF}CD-HIT${NORM}\" using Homebrew" &&
							brew install homebrew/science/cd-hit &&
							echo "\"${REDF}CD-HIT${NORM}\" is available. ${GREF}OK.${NORM}"
								} || {
								echo
								echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of${NORM} \"${REDF}CD-HIT${NORM}\" ${CYAF}failed.${NORM} Please, do it manually. For details see"
								echo "\"${REDF}brew info homebrew/science/cd-hit${NORM}\" and \"${REDF}brew help${NORM}\"."
								echo
								exit 1
								}
							else
								echo "This is not Mac OS X. Going to compile..."
								compilecdhit $SCRIPTDIR/src/cd-hit-v4.6.4-2015-0603
							fi
						break
						;;
					M|m)
						echo
						echo "Please, go to ${REDF}http://weizhongli-lab.org/cd-hit/${NORM}, install \"${REDF}CD-HIT${NORM}\" and ensure"
						echo "  it is in ${BOLD}PATH${NORM}."
						echo
						exit
						;;
					*) echo "${CYAF}Wrong option.${NORM} Use ${REDF}C${NORM}, ${REDF}S${NORM}, ${REDF}B${NORM}, ${REDF}H${NORM} or ${REDF}M${NORM}." && read CDHIT;;
					esac
				done
		else
			exit 1
		fi
	}

# Input files
CHECKFILEREADOUT=""
echo

# Plastome reference in FASTA
readinputfile -c "plastome reference sequence input file in FASTA format" $REFERENCECP0
REFERENCECP0=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Geneious output files are input files here - consensus and unused sequences (TSV)
readinputfile -x "input file in TSV format (output of Geneious assembly)" $TSVLIST
TSVLIST=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Geneious output files are input files here - consensus and unused sequences (FASTA)
readinputfile -z "input file in FASTA format (output of Geneious assembly)" $SEQUENCES0
SEQUENCES0=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input - reference genome - cpDNA
echo "Input file: ${REDF}$REFERENCECP0${NORM}"
# Input - reference genome - cpDNA - checked not to be interleaved - temporary file - will be deleted
REFERENCECP="${REFERENCECP0%.*}_non-interleaved.fasta"
# Input - Geneious output file in TSV
echo "Input file: ${REDF}$TSVLIST${NORM}"
# Input - Geneious output file in FASTA
echo "Input file: ${REDF}$SEQUENCES0${NORM}"
# Input - Geneious output file in FASTA - checked not to be interleaved - temporary file - will be deleted
SEQUENCES="${SEQUENCES0%.*}_non-interleaved.fasta"
# Sequences converted from FASTA to tabular format - temporary file - will be deleted
SEQUENCESTAB="${OUTPUTFILENAME%.*}.tab"
# Assembled sequences in TSV - temporary file - will be deleted
SEQUENCESTABASSE="${OUTPUTFILENAME%.*}_assembled.tab"
# Filtered probes - temporary file - will be deleted
SEQUENCESPROBESLOCUSLENGTH="${OUTPUTFILENAME%.*}_probes_120-600bp.tab"
# List of usable contigs for joining - temporary file - will be deleted
SEQUENCESPROBESLOCUSLENGTHFORJOIN="${OUTPUTFILENAME%.*}_probes_120-600bp_fin_for_join"
# Exons of a certain minimum length - temporary file - will be deleted
SEQUENCESTABASSEBAITL="${OUTPUTFILENAME%.*}_120bp_assembled_less_than_1kb_transcript_fin.tab"
# Sorted exons of a certain minimum length - temporary file - will be deleted
SEQUENCESTABASSEBAITLSORT="${OUTPUTFILENAME%.*}_120bp_assembled_less_than_1kb_transcript_fin_sorted.tab"
# Exons of a certain minimum length and exons making up genes of a certain minimum total locus length - temporary file - will be deleted
SEQUENCESPROBES120600FIN="${OUTPUTFILENAME%.*}_probes_120-600bp_fin.tab"
# Temporal files when converting from TAB to FASTA - temporary file - will be deleted
SEQUENCESPROBES120600MODIF="${OUTPUTFILENAME%.*}_probes_120-600bp_modified_fin.tab"
SEQUENCESPROBES120600ASSEM="${OUTPUTFILENAME%.*}_probes_120-600bp_assembled_fin.fasta"
SEQUENCESPROBES120600CONTIG="${OUTPUTFILENAME%.*}_probes_120-600bp_contig_fin.fasta"
# Preliminary probe sequences - temporary file - will be deleted
PROBEPRELIM0="${OUTPUTFILENAME%.*}_prelim_probe_seq0.fasta"
# Preliminary probe sequences - corrected labels
PROBEPRELIM="${OUTPUTFILENAME%.*}_prelim_probe_seq.fasta"
# Unclustered exons and clustered exons with 100% sequence identity
PROBEPRELIMCLUSTER100="${OUTPUTFILENAME%.*}_prelim_probe_seq_cluster_100.fasta"
# Unclustered exons and clustered exons with more than a certain sequence similarity - temporary file - will be deleted
PROBEPRELIMCLUSTER90="${OUTPUTFILENAME%.*}_prelim_probe_seq_cluster_90.fasta"
# Unclustered exons and clustered exons with more than a certain sequence similarity (CLSTR file)
UNIQUEPROBEPRELIMCLUSTER90="${OUTPUTFILENAME%.*}_unique_prelim_probe_seq_cluster_90.clstr"
# List of unclustered exons / exons with less than a certain sequence similarity (TXT file) - temporary file - will be deleted
UNIQUEPROBEPRELIM="${OUTPUTFILENAME%.*}_unique_prelim_probe_seq.txt"
# Unclustered exons / exons with less than a certain sequence similarity (FASTA file)
UNIQUEPROBEPRELIMF="${OUTPUTFILENAME%.*}_unique_prelim_probe_seq.fasta"
# Sequence similarity checked by CD-HIT - temporary file - will be deleted
PROBEPRELIMCDHIT="${OUTPUTFILENAME%.*}_sim_test.txt"
# Exons of a certain minimum length making up genes of a certain minimum total locus length
PROBEPRELIMCDHIT2="${OUTPUTFILENAME%.*}_similarity_test.fasta"
# Extracted exons making up genes of a certain minimum total locus length - temporary file - will be deleted
PROBEPRELIMFORJOIN="${OUTPUTFILENAME%.*}_similarity_test_assemblies_for_join"
# Modified and sorted exons - temporary file - will be deleted
PROBEPRELIMSORT="${OUTPUTFILENAME%.*}_similarity_test_assemblies_sort.tab"
# Exons of a certain minimum length and exons making up genes of a certain minimum total locus length - temporary file - will be deleted
PROBEPRELIMFIN="${OUTPUTFILENAME%.*}_similarity_test_assemblies_fin.tab"
# Probes in FASTA
PROBESEQUENCES="${OUTPUTFILENAME%.*}_target_enrichment_probe_sequences.fasta"
# Probes in tab - temporary file - will be deleted
PROBESEQUENCESNUM="${OUTPUTFILENAME%.*}_target_enrichment_probe_sequences.tab"
# Putative cpDNA genes in probe set
PROBESEQUENCESCP="${OUTPUTFILENAME%.*}_possible_cp_dna_genes_in_probe_set.pslx"

# Assemble the obtained sequences in contigs

# Check EOL of input files
echo
eolcheck $REFERENCECP0
eolcheck $TSVLIST
eolcheck $SEQUENCES0

# Check if FASTA input files are non-interleaved (required) - if not, FASTA input file converted
echo
echo "Checking if input FASTA files are non-interleaved (required) - interleaved"
echo "  FASTA files are converted not to be interleaved"
echo
noninterleavedfasta $REFERENCECP0 $REFERENCECP
noninterleavedfasta $SEQUENCES0 $SEQUENCES

# Step 8: Retention of those contigs that comprise exons ≥ bait length (default is 120 bp) and have a certain minimum total locus length

echo "${REDF}Step 8 of the pipeline${NORM} - retention of those contigs that comprise exons ≥ bait"
echo "  length (${CYAF}$BAITL${NORM} bp) and have a certain minimum total locus length"

# Check if TSV output of Geneious contains at least required columns
echo
echo "Checking if ${REDF}$TSVLIST${NORM} has all required columns"
REQUIREDCOLS="Required columns in `echo $TSVLIST` are \"`echo ${CYAF}`# Sequences`echo ${NORM}`\",\n  \"`echo ${CYAF}`% Pairwise Identity`echo ${NORM}`\", \"`echo ${CYAF}`Description`echo ${NORM}`\", \"`echo ${CYAF}`Mean Coverage`echo ${NORM}`\", \"`echo ${CYAF}`Name`echo ${NORM}`\"\n  and \"`echo ${CYAF}`Sequence Length`echo ${NORM}`\". Please, export the TSV file again."
echo
if grep -q "# Sequences" $TSVLIST; then
	echo "Column \"${CYAF}# Sequences${NORM}\" is present in ${REDF}$TSVLIST${NORM}. ${GREF}OK.${NORM}"
	if grep -q "Pairwise Identity" $TSVLIST; then
		echo "Column \"${CYAF}Pairwise Identity${NORM}\" is present in ${REDF}$TSVLIST${NORM}. ${GREF}OK.${NORM}"
		if grep -q "Description" $TSVLIST; then
		echo "Column \"${CYAF}Description${NORM}\" is present in ${REDF}$TSVLIST${NORM}. ${GREF}OK.${NORM}"
			if grep -q "Mean Coverage" $TSVLIST; then
			echo "Column \"${CYAF}Mean Coverage${NORM}\" is present in ${REDF}$TSVLIST${NORM}. ${GREF}OK.${NORM}"
				if grep -q "Name" $TSVLIST; then
				echo "Column \"${CYAF}Name${NORM}\" is present in ${REDF}$TSVLIST${NORM}. ${GREF}OK.${NORM}"
					if grep -q "Sequence Length" $TSVLIST; then
					echo "Column \"${CYAF}Sequence Length${NORM}\" is present in ${REDF}$TSVLIST${NORM}. ${GREF}OK.${NORM}"
						else
						echo "${REDF}${BOLD}Error!${NORM} Column \"${CYAF}Sequence Length${NORM}\" is missing!"
						echo -e "$REQUIREDCOLS"
						exit 1
						fi
						else
						echo "${REDF}${BOLD}Error!${NORM} Column \"${CYAF}Name${NORM}\" is missing!"
						echo -e "$REQUIREDCOLS"
						exit 1
						fi
					else
					echo "${REDF}${BOLD}Error!${NORM} Column \"${CYAF}Mean Coverage${NORM}\" is missing!"
					echo -e "$REQUIREDCOLS"
					exit 1
					fi
				else
				echo "${REDF}${BOLD}Error!${NORM} Column \"${CYAF}Description${NORM}\" is missing!"
				echo -e "$REQUIREDCOLS"
				exit 1
				fi
			else
			echo "${REDF}${BOLD}Error!${NORM} Column \"${CYAF}Pairwise Identity${NORM}\" is missing!"
			echo -e "$REQUIREDCOLS"
			exit 1
			fi
		else
		echo "${REDF}${BOLD}Error!${NORM} Column \"${CYAF}# Sequences${NORM}\" is missing!"
		echo -e "$REQUIREDCOLS"
		exit 1
	fi || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking of required columns failed.${NORM} Aborting."
		echo "Check if file ${REDF}$TSVLIST${NORM} is correct TSV file correctly exported from Geneious."
		echo
		exit 1
		}
echo

# Check if TSV output of Geneious contains only required columns or more
if egrep -q "^# Sequences[[:blank:]]+% Pairwise Identity[[:blank:]]+Description[[:blank:]]+Mean Coverage[[:blank:]]+Name[[:blank:]]+Sequence Length$" $TSVLIST
	then
		echo "${REDF}$TSVLIST${NORM} is correct input file. ${GREF}OK.${NORM}"
		TSVLIST2=$TSVLIST
		echo
	else
		echo "Input file ${REDF}$TSVLIST${NORM} seems to contain more columns or columns in"
		echo "  another order than required."
		echo "Needed columns will be extracted in required order."
		$SCRIPTDIR/geneious_column_separator.pl $TSVLIST || {
			echo
			echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Extraction failed.${NORM} Aborting."
			echo "Either script ${REDF}$SCRIPTDIR/geneious_column_separator.pl${NORM}"
			echo "  is missing or there is something wrong with ${REDF}$TSVLIST${NORM}"
			echo "Please, prepare required file manually."
			echo -e "$REQUIREDCOLS"
			echo
			exit 1
			}
		TSVLIST2="${TSVLIST%.*}.columns.tsv"
		echo "File with extracted columns was saved as"
		echo "${REDF}$TSVLIST2${NORM} for possible later usage."
		confirmgo
	fi || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking of required columns failed.${NORM} Aborting."
		echo "Check if file ${REDF}$TSVLIST${NORM} is correct TSV file correctly exported from Geneious."
		echo
		exit 1
		}
echo

# Check the statistics

echo "${REDF}Assembly statistics${NORM}"
echo

# Calculation of the total number of base pairs, based on exons ≥ bait length
echo "${CYAF}Total number of base pairs:${NORM}"
{ cut -f 6 $TSVLIST2 | awk '$1>'"$BAITLN"'' | awk '{s+=$1}END{print s}'; } || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking statistics failed.${NORM} Aborting. Check if file"
	echo "${REDF}$TSVLIST2${NORM} is correct TSV file containing all required columns:"
	echo -e "$REQUIREDCOLS"
	echo
	exit 1
	}
confirmgo

# Check number of contigs
echo "${CYAF}Number of contigs longer than ${REDF}$BAITLN${CYAF}:${NORM}"
{ cut -f 6 $TSVLIST2 | awk '$1>'"$BAITLN"'' | wc -l; } || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking number of contigs failed.${NORM} Aborting. Check if file"
	echo "${REDF}$TSVLIST2${NORM} is correct TSV file containing all required columns"
	echo -e "$REQUIREDCOLS"
	echo
	exit 1
	}
confirmgo

# Convert FASTA to TAB
echo "Converting FASTA to TAB"
fasta2tab $SEQUENCES $SEQUENCESTAB
echo

# Modify labels of FASTA sequences to ensure them to work correctly
{ echo "Checking and modifying FASTA sequence labels" &&
	sed -i 's/:.*\t/\t/' $SEQUENCESTAB &&
	awk -F '[_\t]' '{ printf "%012d_", $1; print; }' $SEQUENCESTAB > $SEQUENCESTAB.temp &&
	mv $SEQUENCESTAB.temp $SEQUENCESTAB &&
	sed -i 's/_[[:digit:]]\+//' $SEQUENCESTAB; } || {
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Modifications of FASTA labels failed.${NORM} Aborting."
		echo "Check if ${REDF}$SEQUENCESTAB${NORM} is correct file."
		echo
		exit 1
		}
echo

# Separate the assembled sequences
echo "Separating assembled sequences"
grep 'Assembly\|Contig' $SEQUENCESTAB > $SEQUENCESTABASSE || {
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Separation of assembled sequences failed.${NORM} Aborting."
	echo "Check if ${REDF}$SEQUENCESTAB${NORM} is correct file (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Retention of those contigs that comprise exons ≥ bait length and have a certain minimum total locus length
# Allowing the values 80, 100, 120 for bait / minimum exon length and 360, 480, 600, 720, 840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920, 2040 for minimum total locus length

echo "${CYAF}Number of assembled sequences:${NORM}"
echo "Length of exons ≥${CYAF}$BAITL${NORM} bp."
echo -e "${REDF}G${CYAF}enes of length${NORM}\t\t${REDF}T${CYAF}otal bp${NORM}\t${REDF}N${CYAF}umber of genes${NORM}"
for LOCUSLENGTH in 0360 0480 0600 0720 0840 0960 1080 1200 1320 1440 1560 1680 1800 1920 2040; do
	LOCUSLENGTHN=$(expr $LOCUSLENGTH - 1)
	echo -e "≥$LOCUSLENGTH bp\t\t$(awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$LOCUSLENGTHN"'' | awk '{s+=$3;c++}END{print s}')\t\t$(awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$LOCUSLENGTHN"'' | wc -l)"
	done || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking number of assembled sequences failed.${NORM} Aborting."
		echo "Check if file ${REDF}$SEQUENCESTABASSE${NORM} is correct (FASTA converted into TAB)."
		echo
		exit 1
		}
echo

# Select the optimal minimum total locus length
if [ "$STARTINI" == "I" ]; then
	echo "${CYAF}Select minimum total locus length.${NORM} Possible values are ${REDF}360${NORM}, ${REDF}480${NORM}, ${REDF}600${NORM},"
	echo "  ${REDF}720${NORM}, ${REDF}840${NORM}, ${REDF}960${NORM}, ${REDF}1080${NORM}, ${REDF}1200${NORM}, ${REDF}1320${NORM}, ${REDF}1440${NORM}, ${REDF}1560${NORM}, ${REDF}1680${NORM}, ${REDF}1800${NORM}, ${REDF}1920${NORM} or ${REDF}2040${NORM}."
	read MINLOCUSLENGTHTEST
	while :
		do
		case "$MINLOCUSLENGTHTEST" in
			360)
				MINLOCUSLENGTH=360
				break
				;;
			480)
				MINLOCUSLENGTH=480
				break
				;;
			600)
				MINLOCUSLENGTH=600
				break
				;;
			720)
				MINLOCUSLENGTH=720
				break
				;;
			840)
				MINLOCUSLENGTH=840
				break
				;;
			960)
				MINLOCUSLENGTH=960
				break
				;;
			1080)
				MINLOCUSLENGTH=1080
				break
				;;
			1200)
				MINLOCUSLENGTH=1200
				break
				;;
			1320)
				MINLOCUSLENGTH=1320
				break
				;;
			1440)
				MINLOCUSLENGTH=1440
				break
				;;
			1560)
				MINLOCUSLENGTH=1560
				break
				;;
			1680)
				MINLOCUSLENGTH=1680
				break
				;;
			1800)
				MINLOCUSLENGTH=1800
				break
				;;
			1920)
				MINLOCUSLENGTH=1920
				break
				;;
			2040)
				MINLOCUSLENGTH=2040
				break
				;;
			*)
				echo "${CYAF}Wrong option.${NORM} Use ${REDF}360${NORM}, ${REDF}480${NORM}, ${REDF}600${NORM}, ${REDF}720${NORM}, ${REDF}840${NORM}, ${REDF}960${NORM}, ${REDF}1080${NORM}, ${REDF}1200${NORM}, ${REDF}1320${NORM}, ${REDF}1440${NORM},"
				echo "  ${REDF}1560${NORM}, ${REDF}1680${NORM}, ${REDF}1800${NORM}, ${REDF}1920${NORM} or ${REDF}2040${NORM}."
				read MINLOCUSLENGTHTEST
				;;
			esac
		done
	fi

echo
echo "${CYAF}Total locus length${NORM} is set to ${REDF}$MINLOCUSLENGTH bp${NORM}."
echo

# Variable to calculate with minimal total locus length
MINLOCUSLENGTHN=$(expr $MINLOCUSLENGTH - 1)

# Saving sequences with selected minimum total locus length
echo "Saving sequences of selected length (≥$MINLOCUSLENGTH bp)"
{ awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' > $SEQUENCESPROBESLOCUSLENGTH; } || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Saving sequences of selected length failed.${NORM} Aborting."
	echo "Check if file ${REDF}$SEQUENCESTABASSE${NORM} is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Create the final FASTA file for the Hyb-Seq probes

# Extract and sort the exons making up genes of a certain minimum total length
echo "Extracting and sorting the exons making up genes of ≥${CYAF}$MINLOCUSLENGTH${NORM} bp"
{ sed 's/^/Assembly_/' $SEQUENCESPROBESLOCUSLENGTH | cut -f 1 -d " " | sort -k 1,1 > $SEQUENCESPROBESLOCUSLENGTHFORJOIN; } || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Extraction and sort of the exons failed.${NORM} Aborting."
	echo "Check if file ${REDF}$SEQUENCESPROBESLOCUSLENGTH${NORM} is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Make a file with all exons of a certain minimum length
echo "Selecting ≥${CYAF}$BAITL${NORM} bp exons"
{ awk '{print $1"\t"length($2)"\t"$2}' $SEQUENCESTABASSE | awk '$2>'"$BAITLN"'' > $SEQUENCESTABASSEBAITL; } || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Selection of the exons failed.${NORM} Aborting."
	echo "Check if file ${REDF}$SEQUENCESTABASSE${NORM} is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Make the assembly number the first field and sort
echo "Sorting exons ≥${CYAF}$BAITL${NORM} bp"
grep '[Cc]ontig' $SEQUENCESTABASSEBAITL | sed 's/^.*\([[:digit:]]\{12\}\).*\([Cc]ontig_[[:digit:]]\{1,\}\).*\>\t\([[:digit:]]\{1,\}\)\t\([[:alpha:]]\{1,\}$\)/Assembly_\1\t\2\t\3\t\4/' | sort -k 1,1 > $SEQUENCESTABASSEBAITLSORT && REMAINING="YES"
# Geneious 9 has different naming scheme of contigs (it does not use keyword "Contig")
if [ ! -s "$SEQUENCESTABASSEBAITLSORT" ]; then
	grep '[Aa]ssembly' $SEQUENCESTABASSEBAITL | sed 's/^.*\([[:digit:]]\{12\}\).*[Aa]ssembly_\([[:digit:]]\{1,\}\).*\>\t\([[:digit:]]\{1,\}\)\t\([[:alpha:]]\{1,\}$\)/Assembly_\1\t\2\t\3\t\4/' | sort -k 1,1 > $SEQUENCESTABASSEBAITLSORT
	REMAINING="NO"
	fi
echo

# Make a file with all exons of a certain minimum length and making up genes of a certain minimum length
echo "Selecting all exons ≥${CYAF}$BAITL${NORM} bp and all exons making up genes of ≥${CYAF}$MINLOCUSLENGTH${NORM} bp"
join $SEQUENCESPROBESLOCUSLENGTHFORJOIN $SEQUENCESTABASSEBAITLSORT > $SEQUENCESPROBES120600FIN || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Selection of the exons failed.${NORM} Aborting."
	echo "Check if files ${REDF}$SEQUENCESPROBESLOCUSLENGTHFORJOIN${NORM} and ${REDF}$SEQUENCESTABASSEBAITLSORT${NORM}"
	echo "  are correct (tabular files listing respective exons and their sequences)."
	echo
	exit 1
	}
echo

# Convert TAB to FASTA
echo "Converting TAB to FASTA"
{ sed 's/ /_/' $SEQUENCESPROBES120600FIN | sed 's/ /_/' > $SEQUENCESPROBES120600MODIF &&
	sed 's/^/>/' $SEQUENCESPROBES120600MODIF | sed 's/ /\n/' > $SEQUENCESPROBES120600ASSEM; } || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of FASTA to TAB failed.${NORM} Aborting."
		echo "Check if file ${REDF}$SEQUENCESPROBES120600FIN${NORM} is correct FASTA file."
		echo
		exit 1
		}
echo

# Remaining assemblies have to be selected and added to the FASTA file of the probes
# NOTE: Geneious 9 does not use keyword "Contig"
if [ "$REMAINING"=="YES" ]; then
	grep -v '[Cc]ontig' $SEQUENCESTABASSEBAITL | awk '$2>'"$MINLOCUSLENGTHN"'' | sed 's/^/>/' | sed 's/\t/_/' | sed 's/\t/\n/' > $SEQUENCESPROBES120600CONTIG
	fi
echo

# Combine the two FASTA files
echo "Writing FASTA file with preliminary probe sequences"
cat $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG > $PROBEPRELIM0 || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Writing of preliminary probe sequences failed.${NORM} Aborting."
	echo "Check if files ${REDF}$SEQUENCESPROBES120600ASSEM${NORM} and ${REDF}$SEQUENCESPROBES120600CONTIG${NORM} are correct FASTA files."
	echo
	exit 1
	}
echo

# Ensure all sequences have correct labels
echo "Ensuring all sequences have correct labels"
sed 's/^>[^0123456789]*\([[:digit:]]\{12\}\)[^0123456789]*\([[:digit:]]\{1,\}\)[^0123456789]*\([[:digit:]]\{1,\}\)$/>Assembly_\1_Contig_\2_\3/' $PROBEPRELIM0 > $PROBEPRELIM || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking of FASTA sequence labels failed.${NORM} Aborting."
	echo "Check if file ${REDF}$PROBEPRELIM0${NORM} is correct FASTA file."
	echo
	exit 1
	}
echo
echo "${CYAF}Preliminary probe sequences${NORM} saved as"
echo "  ${REDF}$PROBEPRELIM${NORM} for possible later usage."
confirmgo

# Step 9: Make the final quality control of the probe sequences

echo "${REDF}Step 9 of the pipeline${NORM} - removal of probe sequences sharing ≥90% sequence"
echo "  similarity"
echo

# Check for sequence similarity between the developed probe sequences with CD-HIT-EST

# Clustering exons with 100% sequence identity: retaining unclustered exons and, in case of 100% sequence identity, retaining the longest exon
echo "${CYAF}Checking sequence similarity between the probe sequences (exons)${NORM}"
echo "  Detecting identical probe sequences and retaining the longest one in such a"
echo "  case. Retaining also the unclustered probe sequences."
echo
cd-hit-est -i $PROBEPRELIM -o $PROBEPRELIMCLUSTER100 -d 0 -c 1.0 -p 1 || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking of the sequence similarity failed.${NORM} Aborting."
	echo "Check if file ${REDF}$PROBEPRELIM${NORM} is correct FASTA file."
	echo
	exit 1
	}
echo
echo "Clustered exons with 100% sequence identity were saved as"
echo "${REDF}$PROBEPRELIMCLUSTER100${NORM} for possible later usage."
confirmgo

# Clustering and removing exons with more than a certain sequence similarity
echo "${CYAF}Detecting and removing probe sequences (exons)${NORM} that are similar to each other"
echo "  above a certain threshold"
echo
cd-hit-est -i $PROBEPRELIMCLUSTER100 -o $PROBEPRELIMCLUSTER90 -d 0 -c $CDHITSIM -p 1 -g 1 || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking of the probe sequence failed.${NORM} Aborting."
	echo "Check if file ${REDF}$PROBEPRELIMCLUSTER100${NORM} is correct FASTA file."
	echo
	exit 1
	}
echo
# Finding those clusters from a CD-HIT CLSTR file that include only one sequence or multiple sequences with 100% identity (in which case the longest sequence is choosen)
python $SCRIPTDIR/grab_singleton_clusters.py -i $PROBEPRELIMCLUSTER90.clstr -o $UNIQUEPROBEPRELIMCLUSTER90 || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Checking of the probe sequence failed.${NORM} Aborting."
	echo "Check if file ${REDF}$PROBEPRELIMCLUSTER90.clstr${NORM} is correct output of ${CYAF}cd-hit-est${NORM}."
	echo
	exit 1
	}
echo
echo "Unclustered exons and clustered exons with 100% identity were saved as"
echo "  ${REDF}$UNIQUEPROBEPRELIMCLUSTER90${NORM} for possible later usage."
confirmgo

echo "Postprocessing extracted sequences"
{ grep -v '>Cluster' $UNIQUEPROBEPRELIMCLUSTER90 | cut -d ' ' -f 2 | sed -e 's/\.\.\./\\\>/' -e 's/^/^/' > $UNIQUEPROBEPRELIM &&
grep -A 1 -f $UNIQUEPROBEPRELIM $PROBEPRELIMCLUSTER100 | sed '/^--$/d' > $UNIQUEPROBEPRELIMF; } || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Postprocessing of extracted sequences failed.${NORM} Aborting."
	echo "Check if file ${REDF}$UNIQUEPROBEPRELIMCLUSTER90${NORM} is correct output of ${CYAF}cd-hit-est${NORM}."
	echo
	exit 1
	}
echo
echo "Postprocessed extracted sequences were saved as"
echo "  ${REDF}$UNIQUEPROBEPRELIMF${NORM} for possible later usage."
confirmgo

# Step 10: Retention of those probe sequences that comprise exons of a certain minimum length (default is 120 bp) and have a certain minimum total locus length

echo "${REDF}Step 10 of the pipeline${NORM} - retention of those probe sequences that comprise exons"
echo "  of a certain minimum length (${CYAF}$BAITL${NORM} bp) and have a certain minimum total locus length"
echo

# One of the three outfiles of last steps of part 9 is a FASTA file, it has to be converted to TAB
echo "Converting FASTA from step 9 to TAB"
fasta2tab $UNIQUEPROBEPRELIMF $PROBEPRELIMCDHIT
echo

# # Count the exons of a certain minimum length (default ≥120 bp)
# echo "${CYAF}Number of exons of a minimum length ≥$BAITL bp:${NORM}"
# awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT | awk '{s+=$2;c++}END{print s}'
# confirmgo
#
# # Count the exons of a certain minimum length making up genes of a certain minimum length
# echo "$(awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT | sed 's/_/\t/g' | cut -f 2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' | awk '{s+=$3;c++}END{print s}') ${CYAF}of the assemblies making up genes of ≥$MINLOCUSLENGTH bp,"
# echo "  comprised of ${REDF}$(awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT | sed 's/_/\t/g' | cut -f 2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' | wc -l)${NORM} ${CYAF}putative exons ≥${REDF}$BAITL${NORM} ${CYAF}bp${NORM}."
# confirmgo

echo "Writing the exons into temporal file"
awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT | sed 's/_/\t/g' | cut -f 2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' > $PROBEPRELIMCDHIT2
echo

# Extract and sort the exons making up genes of a minimum total length
echo "Extracting and sorting the exons making up genes of ≥$MINLOCUSLENGTH bp"
sed 's/^/Assembly_/' $PROBEPRELIMCDHIT2 | cut -f 1 -d " " | sort -k 1,1 > $PROBEPRELIMFORJOIN
echo

# Modify the exon number and sort
echo "Modifying the exon number and sorting"
sed 's/_C/\tC/' $PROBEPRELIMCDHIT | sort -k 1,1 > $PROBEPRELIMSORT
echo

# Make a file with all exons of a certain minimum length making up genes of a certain minimum total length
echo "Joining all exons ≥${CYAF}$BAITL${NORM} bp and making up genes of ≥$MINLOCUSLENGTH bp"
join $PROBEPRELIMFORJOIN $PROBEPRELIMSORT > $PROBEPRELIMFIN
echo

# # Convert TAB to FASTA
# echo "Converting TAB to FASTA"
# sed 's/^\(.\+\) \(Contig\)/>\1_\2/' $PROBEPRELIMFIN | sed 's/ /\n/' > $PROBESEQUENCES
# echo

# Probe design summary

# Calculation of the total number of base pairs

echo "Calculating the total number of base pairs"
# echo "Converting FASTA to TAB"
# fasta2tab $PROBESEQUENCES $PROBESEQUENCESNUM

echo "$(awk '{print $1"\t"length($2)}' $PROBEPRELIMFIN | sed 's/_/\t/g' | cut -f 2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' | awk '{s+=$3;c++}END{print s}') ${CYAF}bp make up genes of ≥$MINLOCUSLENGTH bp,"

# echo
# echo "${CYAF}Total number of base pairs:${NORM} $(awk '{print $1"\t"length($2)}' $PROBESEQUENCESNUM | awk '{s+=$2;c++}END{print s}')"
# echo

# Calculation of the total number of genes

echo "Calculating the total number of genes"
echo " There are ${REDF}$(awk '{print $1"\t"length($2)}' $PROBEPRELIMFIN | sed 's/_/\t/g' | cut -f 2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' | wc -l)${NORM} ${CYAF}genes in total."
confirmgo

# Calculation of the total number of exons

echo "Calculating the total number of exons"
echo "${CYAF}Total number of exons:${NORM} $(wc -l $PROBEPRELIMFIN) ≥${REDF}$BAITL${NORM} ${CYAF}bp${NORM}."
confirmgo

# Convert TAB to FASTA
echo "Converting TAB to FASTA"
sed 's/^\(.\+\) \(Contig\)/>\1_\2/' $PROBEPRELIMFIN | sed 's/ /\n/' > $PROBESEQUENCES
echo

# # Calculation of the total number of genes
# echo "Calculating the total number of genes"
# echo
# echo -e "${REDF}G${CYAF}enes of length${NORM}:\t≥$LOCUSLENGTH bp"
# echo -e "${REDF}T${CYAF}otal bp${NORM}:\t\t$(awk '{print $1"\t"length($2)}' $PROBESEQUENCESNUM | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$LOCUSLENGTHN"'' | awk '{s+=$3;c++}END{print s}')"
# echo -e "${REDF}N${CYAF}umber of genes${NORM}:\t$(awk '{print $1"\t"length($2)}' $PROBESEQUENCESNUM | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$LOCUSLENGTHN"'' | wc -l)"
# confirmgo

echo "${REDF}${BOLD}Success!${NORM}"
echo
echo "${BLUF}================================================================================${NORM}"
echo "${CYAF}Final output file was written as${NORM}"
echo "${REDF}${BOLD}$PROBESEQUENCES${NORM}"
echo "${CYAF}This file contains the probe sequences.${NORM}"
echo "${BLUF}================================================================================${NORM}"
confirmgo

# Step 11 - removal of putative cpDNA sequences in final probe list

echo "${REDF}Step 11 of the pipeline${NORM} - detection of probe sequences sharing ${CYAF}≥$BLATIDENT%${NORM} sequence"
echo "similarity with the plastome reference"
echo

# Remove remaining cp genes from probe set
echo "${CYAF}Detecting remaining plastid genes in probe set${NORM}"
blat -t=dna -q=dna -minIdentity=$BLATIDENT -out=pslx $REFERENCECP $PROBESEQUENCES $PROBESEQUENCESCP
echo

echo "${BLUF}================================================================================${NORM}"
echo "File ${REDF}${BOLD}$PROBESEQUENCESCP${NORM}"
echo "contains ${CYAF}putative plastid genes in final probe set${NORM}."
echo "We ${REDF}STRONGLY RECOMMEND${NORM} to ${CYAF}remove those genes${NORM} from the final probe set in file"
echo "${REDF}${BOLD}$PROBESEQUENCES${NORM}."
echo "${BLUF}================================================================================${NORM}"
confirmgo

# Remove temporal files
echo "Removing unneeded temporal files"
rm $REFERENCECP $SEQUENCES $SEQUENCESTAB $SEQUENCESTABASSE $SEQUENCESPROBESLOCUSLENGTH $SEQUENCESPROBESLOCUSLENGTHFORJOIN $SEQUENCESTABASSEBAITL $SEQUENCESTABASSEBAITLSORT $SEQUENCESPROBES120600FIN $SEQUENCESPROBES120600MODIF $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG $PROBEPRELIM0 $PROBEPRELIMCLUSTER90 $UNIQUEPROBEPRELIM $PROBEPRELIMCLUSTER100 $UNIQUEPROBEPRELIMF $PROBEPRELIMCDHIT $PROBEPRELIMFORJOIN $PROBEPRELIMSORT $PROBEPRELIMFIN  || { # $PROBESEQUENCESNUM currently not in use
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removal of temporal files failed.${NORM} Remove following files manually:"
	echo "  \"$REFERENCECP\", \"$SEQUENCES\", \"$SEQUENCESTAB\","
	echo "  \"$SEQUENCESTABASSE\", \"$SEQUENCESPROBESLOCUSLENGTH\",\"$SEQUENCESPROBESLOCUSLENGTHFORJOIN\","
	echo "  \"$SEQUENCESTABASSEBAITL\", \"$SEQUENCESTABASSEBAITLSORT\",\"$SEQUENCESPROBES120600FIN\","
	echo "  \"$SEQUENCESPROBES120600MODIF\", \"$SEQUENCESPROBES120600ASSEM\", \"$SEQUENCESPROBES120600CONTIG\","
	echo "  \"$PROBEPRELIM0\", \"$PROBEPRELIMCLUSTER90\", \"$UNIQUEPROBEPRELIM\","
	echo "  \"$PROBEPRELIMCDHIT\", \"$PROBEPRELIMFORJOIN\", \"$PROBEPRELIMSORT\","
	echo "  \"$PROBEPRELIMFIN\" and \"$PROBESEQUENCESNUM\"."
	confirmgo
	}

# List kept files which user can use for another analysis
echo
echo "${CYAF}Following files are kept for possible later usage (see manual for details):${NORM}"
echo "${BLUF}================================================================================${NORM}"
echo "${CYAF}1)${NORM} Preliminary probe sequences:"
echo "${REDF}$PROBEPRELIM${NORM}"
echo "${CYAF}2)${NORM} Unclustered exons and clustered exons with 100% sequence identity:"
echo "${REDF}$PROBEPRELIMCLUSTER100${NORM}"
echo "${CYAF}3)${NORM} Clustered exons with 100% sequence identity:"
echo "${REDF}$UNIQUEPROBEPRELIMCLUSTER90${NORM}"
echo "${CYAF}4)${NORM} Unclustered exons / exons with less than a certain sequence similarity:"
echo "${REDF}$UNIQUEPROBEPRELIMF${NORM}"
echo "${CYAF}5)${NORM} Contigs comprising exons ≥ bait length having a certain minimum total locus length:"
echo "${REDF}$PROBEPRELIMCDHIT2${NORM}"
echo "${CYAF}6)${NORM} Putative plastid genes in final probe set:"
echo "${REDF}$PROBESEQUENCESCP${NORM}"
echo "${CYAF}7)${NORM} Final probe sequences in FASTA format:"
echo "${REDF}$PROBESEQUENCES${NORM}"
echo "${BLUF}================================================================================${NORM}"
confirmgo

echo "Script exited successfully..."
echo

exit

 || {
	echo
	echo "${REDF}${BOLD}Error!${NORM} ${CYAF}XXX failed.${NORM} Aborting."
	echo "Check if file ${REDF}$XXX${NORM} is correct XXX file."
	echo
	exit 1
	}