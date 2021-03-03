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
	echo "Fatal error! Unable to load file \"sondovac_functions\" with required functions! It must be in same directory as \"$0\" Check it and, if needed, download again whole script from https://github.com/V-Z/sondovac/"
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
# Minimum bait/exon length
BAITL=120
# CD-HIT sequence similarity
CDHITSIM=0.9
# BLAT -minIdentity between the probe sequences and the plastome reference
BLATIDENT=90
# Minimum total locus length
MINLOCUSLENGTH=600
# Default name of output files
OUTPUTFILENAME=`realpath "output"`
# Add also unassembled sequences longer than minimal length
REMAINING="YES"

# Create empty variables for file names
REFERENCECP=""
REFERENCECP0=""
TSVLIST=""
SEQUENCES=""
SEQUENCES0=""

# Parse initial arguments
while getopts "hvlrpeo:inc:x:z:b:d:y:k:" START; do
	case "$START" in
		h|v)
			generaloptions
			echo
			echo -e "\tIf options -c, -x and/or -z are used and script is running in interactive mode, those values will be used as defaults, but may be later overwritten."
			echo
			echo -e "\tOptions required for running in non-interactive mode:"
			echo -e "\t-c\tPlastome reference sequence input file in FASTA format."
			echo -e "\t-x\tInput file in TSV format (output of Geneious assembly)."
			echo -e "\t-z\tInput file in FASTA format (output of Geneious assembly)."
			echo
			echo -e "\tOther optional arguments (if not provided, default values are used):"
			echo -e "\t-b\tBait length. Default value: 120 (preferred length for phylogeny, use any of values 80, 100 or 120)."
			echo -e "\t-d\tSequence similarity between the developed probe sequences (parameter \"-c\" of cd-hit-est, see its manual for details). Default value: 0.9 (use decimal number ranging from 0.85 to 0.95)."
			echo -e "\t-y\tSequence similarity between the probes and plastome reference searching for possible plastid genes in probe set (parameter \"-minIdentity\" of BLAT, see its manual for details). Default value: 90 (integer ranging from 85 to 95)."
			echo -e "\t-k\tMinimum total locus length. Default value: 600. Allowed values are 360, 480, 600, 720, 840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920 and 2040. When running in interactive mode, the user will be asked which value to use. A table summarizing the total number of LCN loci and the total number of base pairs for these values will be displayed to facilitate this choice."
			echo -e "\tWARNING! If parameters -b, -d or -y are not provided, default values are taken, and it is not possible to change them later (not even in interactive mode)."
			echo
			exit 2
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
			OUTPUTFILENAME=`realpath $OPTARG`
			echo "Output files will start name with $OUTPUTFILENAME"
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
			REFERENCECP0=$OPTARG
			echo "Plastome reference: $REFERENCECP"
			;;
		x)
			TSVLIST=$OPTARG
			echo "Input file: $TSVLIST"
			;;
		z)
			SEQUENCES0=$OPTARG
			echo "Input file: $SEQUENCES0"
			;;
		b)
			BAITL=$OPTARG
			# Check if provided value makes sense
			case "$BAITL" in
				80) BAITL=80;;
				100) BAITL=100;;
				120) BAITL=120;;
				*) echo
					echo "Error! For parameter \"-b\" you did not provide any of values 80, 100 or 120!"
					echo
					exit 1
				esac
			echo "Bait length: $BAITL"
			;;
		d)
			CDHITSIM=$OPTARG
			# Check if provided value makes sense
			if [ "$(echo 0.85 '<=' $CDHITSIM | bc -l)" = 1 ] && [ "$(echo $CDHITSIM '<=' 0.95 | bc -l)" = 1 ]; then
				echo "Sequence similarity: $CDHITSIM"
				else
					echo
					echo "Error! For parameter \"-d\" you did not provide decimal number ranging from 0.85 to 0.95!"
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
					echo "Error! For parameter \"-y\" you did not provide an integer of range from 85 to 95!"
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
					echo "Error! For parameter \"-k\" you did not provide any of values 360, 480, 600, 720, 840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920 or 2040!"
					echo
					exit 1
				esac
			echo "Minimum exon length: $MINLOCUSLENGTH"
			;;
		?)
			echo
			echo "Invalid option(s)! See \"$0 -h\" for usage options."
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

# Set variables for working directory and PATH
workdirpath

# Check availability of all needed binaries

# Check if realpath is available
checktools realpath

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

# Check if cd-hit-est is available
{ command -v cd-hit-est >/dev/null 2>&1 && echo "\"cd-hit-est\" is available. OK."; } || {
	echo "CD-HIT (command cd-hit-est) is required but not installed or available in PATH. See https://github.com/weizhongli/cdhit for installation of CD-HIT."
	}

# Input files
CHECKFILEREADOUT=""

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
echo

# Input - reference genome - cpDNA
echo "Input file: $REFERENCECP0"
# Input - reference genome - cpDNA - checked not to be interleaved - temporary file - will be deleted
REFERENCECP="${REFERENCECP0%.*}_non-interleaved.fasta"
# Input - Geneious output file in TSV
echo "Input file: $TSVLIST"
# Input - Geneious output file in FASTA
echo "Input file: $SEQUENCES0"
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
# Probes in FASTA (with putative cpDNA genes)
PROBESEQUENCES="${OUTPUTFILENAME%.*}_target_enrichment_probe_sequences_with_pt.fasta"
# Putative cpDNA genes in probe set
PROBESEQUENCESCP="${OUTPUTFILENAME%.*}_possible_cp_dna_genes_in_probe_set.pslx"
# List of names putative cpDNA genes in probe set
PROBESEQUENCESCPLIST="${OUTPUTFILENAME%.*}_possible_cp_dna_genes_in_probe_set.txt"
# Final probes in FASTA (without putative cpDNA genes)
PROBESEQUENCESNOCP="${OUTPUTFILENAME%.*}_target_enrichment_probe_sequences.fasta"

# Assemble the obtained sequences in contigs

# Check EOL of input files
echo
eolcheck $REFERENCECP0
eolcheck $TSVLIST
eolcheck $SEQUENCES0

# Check if FASTA input files are non-interleaved (required) - if not, FASTA input file converted
echo "Checking if input FASTA files are non-interleaved (required) - interleaved FASTA files are converted not to be interleaved"
echo
noninterleavedfasta $REFERENCECP0 $REFERENCECP
noninterleavedfasta $SEQUENCES0 $SEQUENCES

# Step 8: Retention of those contigs that comprise exons ≥ bait length (default is 120 bp) and have a certain minimum total locus length

echo "Step 8 of the pipeline - retention of those contigs that comprise exons ≥ bait length ($BAITL bp) and have a certain minimum total locus length"

# Check if TSV output of Geneious contains at least required columns
echo
echo "Checking if $TSVLIST has all required columns"
REQUIREDCOLS="Required columns in `echo $TSVLIST` are \"`echo `# Sequences`echo `\",\n  \"`echo `% Pairwise Identity`echo `\", \"`echo `Description`echo `\", \"`echo `Mean Coverage`echo `\", \"`echo `Name`echo `\"\n  and \"`echo `Sequence Length`echo `\". Please, export the TSV file again."
echo
if grep -q "# Sequences" $TSVLIST; then
	echo "Column \"# Sequences\" is present in $TSVLIST. OK."
	if grep -q "Pairwise Identity" $TSVLIST; then
		echo "Column \"Pairwise Identity\" is present in $TSVLIST. OK."
		if grep -q "Description" $TSVLIST; then
		echo "Column \"Description\" is present in $TSVLIST. OK."
			if grep -q "Mean Coverage" $TSVLIST; then
			echo "Column \"Mean Coverage\" is present in $TSVLIST. OK."
				if grep -q "Name" $TSVLIST; then
				echo "Column \"Name\" is present in $TSVLIST. OK."
					if grep -q "Sequence Length" $TSVLIST; then
					echo "Column \"Sequence Length\" is present in $TSVLIST. OK."
						else
						echo "Error! Column \"Sequence Length\" is missing!"
						echo -e "$REQUIREDCOLS"
						exit 1
						fi
						else
						echo "Error! Column \"Name\" is missing!"
						echo -e "$REQUIREDCOLS"
						exit 1
						fi
					else
					echo "Error! Column \"Mean Coverage\" is missing!"
					echo -e "$REQUIREDCOLS"
					exit 1
					fi
				else
				echo "Error! Column \"Description\" is missing!"
				echo -e "$REQUIREDCOLS"
				exit 1
				fi
			else
			echo "Error! Column \"Pairwise Identity\" is missing!"
			echo -e "$REQUIREDCOLS"
			exit 1
			fi
		else
		echo "Error! Column \"# Sequences\" is missing!"
		echo -e "$REQUIREDCOLS"
		exit 1
	fi || {
		echo
		echo "Error! Checking of required columns failed. Aborting. Check if file $TSVLIST is correct TSV file correctly exported from Geneious."
		echo
		exit 1
		}
echo

# Check if TSV output of Geneious contains only required columns or more
if egrep -q "^# Sequences[[:blank:]]+% Pairwise Identity[[:blank:]]+Description[[:blank:]]+Mean Coverage[[:blank:]]+Name[[:blank:]]+Sequence Length$" $TSVLIST
	then
		echo "$TSVLIST is correct input file. OK."
		TSVLIST2=$TSVLIST
		echo
	else
		echo "Input file $TSVLIST seems to contain more columns or columns in another order than required. Needed columns will be extracted in required order."
		checktools perl
		$SCRIPTDIR/geneious_column_separator.pl $TSVLIST || {
			echo
			echo "Error! Extraction failed. Aborting. Either script $SCRIPTDIR/geneious_column_separator.pl is missing or there is something wrong with $TSVLIST Please, prepare required file manually."
			echo -e "$REQUIREDCOLS"
			echo
			exit 1
			}
		TSVLIST2="${TSVLIST%.*}.columns.tsv"
		echo "File with extracted columns was saved as $TSVLIST2 for possible later usage."
		confirmgo
	fi || {
		echo
		echo "Error! Checking of required columns failed. Aborting. Check if file $TSVLIST is correct TSV file correctly exported from Geneious."
		echo
		exit 1
		}

# Check the statistics

echo "Assembly statistics"
echo

# Calculation of the total number of base pairs, based on exons ≥ bait length
{ echo "Total number of base pairs: `cut -f 6 $TSVLIST2 | awk '$1>'"$BAITLN"'' | awk '{s+=$1}END{print s}'`."; } || {
	echo
	echo "Error! Checking statistics failed. Aborting. Check if file $TSVLIST2 is correct TSV file containing all required columns:"
	echo -e "$REQUIREDCOLS"
	echo
	exit 1
	}
confirmgo

# Check number of contigs
{ echo "Number of contigs longer than $BAITL bp: `cut -f 6 $TSVLIST2 | awk '$1>'"$BAITLN"'' | wc -l`."; } || {
	echo
	echo "Error! Checking number of contigs failed. Aborting. Check if file $TSVLIST2 is correct TSV file containing all required columns"
	echo -e "$REQUIREDCOLS"
	echo
	exit 1
	}
confirmgo

# Convert FASTA to TAB
echo "Converting FASTA to TAB"
fasta2tab $SEQUENCES $SEQUENCESTAB

# Modify labels of FASTA sequences to ensure them to work correctly
{ echo "Checking and modifying FASTA sequence labels" &&
	sed -i 's/:.*\t/\t/' $SEQUENCESTAB &&
	awk -F '[_\t]' '{ printf "%012d_", $1; print; }' $SEQUENCESTAB > $SEQUENCESTAB.temp &&
	mv $SEQUENCESTAB.temp $SEQUENCESTAB &&
	sed -i 's/_[[:digit:]]\+//' $SEQUENCESTAB; } || {
		echo "Error! Modifications of FASTA labels failed. Aborting. Check if $SEQUENCESTAB is correct file."
		echo
		exit 1
		}
echo

# Separate the assembled sequences
echo "Separating assembled sequences"
grep 'Assembly\|Contig' $SEQUENCESTAB > $SEQUENCESTABASSE || {
	echo "Error! Separation of assembled sequences failed. Aborting. Check if $SEQUENCESTAB is correct file (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Retention of those contigs that comprise exons ≥ bait length and have a certain minimum total locus length
# Allowing the values 80, 100, 120 for bait / minimum exon length and 360, 480, 600, 720, 840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920, 2040 for minimum total locus length

echo "Number of assembled sequences:"
echo "Length of exons ≥$BAITL bp."
echo -e "Genes of length\t\tTotal bp\tNumber of genes"
for LOCUSLENGTH in 0360 0480 0600 0720 0840 0960 1080 1200 1320 1440 1560 1680 1800 1920 2040; do
	LOCUSLENGTHN=$(expr $LOCUSLENGTH - 1)
	echo -e "≥$LOCUSLENGTH bp\t\t$(awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$LOCUSLENGTHN"'' | awk '{s+=$3;c++}END{print s}')\t\t$(awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$LOCUSLENGTHN"'' | wc -l)"
	done || {
		echo
		echo "Error! Checking number of assembled sequences failed. Aborting. Check if file $SEQUENCESTABASSE is correct (FASTA converted into TAB)."
		echo
		exit 1
		}
echo

# Select the optimal minimum total locus length
if [ "$STARTINI" == "I" ]; then
	echo "Select minimum total locus length. Possible values are 360, 480, 600, 720, 840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920 or 2040."
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
				echo "Wrong option. Use 360, 480, 600, 720, 840, 960, 1080, 1200, 1320, 1440, 1560, 1680, 1800, 1920 or 2040."
				read MINLOCUSLENGTHTEST
				;;
			esac
		done
	fi

echo "Total locus length is set to $MINLOCUSLENGTH bp."
echo

# Variable to calculate with minimal total locus length
MINLOCUSLENGTHN=$(expr $MINLOCUSLENGTH - 1)

# Saving sequences with selected minimum total locus length
echo "Saving sequences of selected length (≥$MINLOCUSLENGTH bp)"
{ awk '{print $1"\t"length($2)}' $SEQUENCESTABASSE | sed 's/^.*\([[:digit:]]\{12\}\).*\t/\1\t/' | awk '$2>'"$BAITLN"'' | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' > $SEQUENCESPROBESLOCUSLENGTH; } || {
	echo
	echo "Error! Saving sequences of selected length failed. Aborting. Check if file $SEQUENCESTABASSE is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Create the final FASTA file for the Hyb-Seq probes

# Extract and sort the exons making up genes of a certain minimum total length
echo "Extracting and sorting the exons making up genes of ≥$MINLOCUSLENGTH bp"
{ sed 's/^/Assembly_/' $SEQUENCESPROBESLOCUSLENGTH | cut -f 1 -d " " | sort -k 1,1 > $SEQUENCESPROBESLOCUSLENGTHFORJOIN; } || {
	echo
	echo "Error! Extraction and sort of the exons failed. Aborting. Check if file $SEQUENCESPROBESLOCUSLENGTH is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Make a file with all exons of a certain minimum length
echo "Selecting ≥$BAITL bp exons"
{ awk '{print $1"\t"length($2)"\t"$2}' $SEQUENCESTABASSE | awk '$2>'"$BAITLN"'' > $SEQUENCESTABASSEBAITL; } || {
	echo
	echo "Error! Selection of the exons failed. Aborting. Check if file $SEQUENCESTABASSE is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
echo

# Make the assembly number the first field and sort
echo "Sorting exons ≥$BAITL bp"
{ grep '[Cc]ontig' $SEQUENCESTABASSEBAITL | sed 's/^.*\([[:digit:]]\{12\}\).*\([Cc]ontig_[[:digit:]]\{1,\}\).*\>\t\([[:digit:]]\{1,\}\)\t\([[:alpha:]]\{1,\}$\)/Assembly_\1\t\2\t\3\t\4/' | sort -k 1,1 > $SEQUENCESTABASSEBAITLSORT && REMAINING="YES"; } || {
	echo
	echo "Error! Sorting of exons failed. Aborting. Check if file $SEQUENCESTABASSEBAITL is correct (FASTA converted into TAB)."
	echo
	exit 1
	}
# Geneious 9 has different naming scheme of contigs (it does not use keyword "Contig")
if [ ! -s "$SEQUENCESTABASSEBAITLSORT" ]; then
	{ grep '[Aa]ssembly' $SEQUENCESTABASSEBAITL | sed 's/^.*\([[:digit:]]\{12\}\).*[Aa]ssembly_\([[:digit:]]\{1,\}\).*\>\t\([[:digit:]]\{1,\}\)\t\([[:alpha:]]\{1,\}$\)/Assembly_\1\t\2\t\3\t\4/' | sort -k 1,1 > $SEQUENCESTABASSEBAITLSORT && REMAINING="NO"; } || {
		echo
		echo "Error! Sorting of exons failed. Aborting. Check if file $SEQUENCESTABASSEBAITL is correct (FASTA converted into TAB)."
		echo
		exit 1
		}
	fi
echo

# Make a file with all exons of a certain minimum length and making up genes of a certain minimum length
echo "Selecting all exons ≥$BAITL bp and all exons making up genes of ≥$MINLOCUSLENGTH bp"
join $SEQUENCESPROBESLOCUSLENGTHFORJOIN $SEQUENCESTABASSEBAITLSORT > $SEQUENCESPROBES120600FIN || {
	echo
	echo "Error! Selection of the exons failed. Aborting. Check if files $SEQUENCESPROBESLOCUSLENGTHFORJOIN and $SEQUENCESTABASSEBAITLSORT are correct (tabular files listing respective exons and their sequences)."
	echo
	exit 1
	}
echo

# Convert TAB to FASTA
echo "Converting TAB to FASTA"
{ sed 's/ /_/' $SEQUENCESPROBES120600FIN | sed 's/ /_/' > $SEQUENCESPROBES120600MODIF &&
	sed 's/^/>/' $SEQUENCESPROBES120600MODIF | sed 's/ /\n/' > $SEQUENCESPROBES120600ASSEM; } || {
		echo
		echo "Error! Conversion of FASTA to TAB failed. Aborting. Check if file $SEQUENCESPROBES120600FIN is correct FASTA file."
		echo
		exit 1
		}
echo

# Remaining assemblies have to be selected and added to the FASTA file of the probes
# Geneious 9 does not use keyword "Contig"
if [ "$REMAINING"=="YES" ]; then
	{ grep -v '[Cc]ontig' $SEQUENCESTABASSEBAITL | awk '$2>'"$MINLOCUSLENGTHN"'' | sed 's/^/>/' | sed 's/\t/_/' | sed 's/\t/\n/' > $SEQUENCESPROBES120600CONTIG; } || {
		echo
		echo "Error! Extraction of remaining exons failed. Aborting. Check if file $SEQUENCESTABASSEBAITL is correct (FASTA converted into TAB)."
		echo
		exit 1
		}
	fi

# Combine the two FASTA files
echo "Writing FASTA file with preliminary probe sequences"
cat $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG > $PROBEPRELIM0 || {
	echo
	echo "Error! Writing of preliminary probe sequences failed. Aborting. Check if files $SEQUENCESPROBES120600ASSEM and $SEQUENCESPROBES120600CONTIG are correct FASTA files."
	echo
	exit 1
	}
echo

# Ensure all sequences have correct labels
echo "Ensuring all sequences have correct labels"
sed 's/^>[^0123456789]*\([[:digit:]]\{12\}\)[^0123456789]*\([[:digit:]]\{1,\}\)[^0123456789]*\([[:digit:]]\{1,\}\)$/>Assembly_\1_Contig_\2_\3/' $PROBEPRELIM0 > $PROBEPRELIM || {
	echo
	echo "Error! Checking of FASTA sequence labels failed. Aborting. Check if file $PROBEPRELIM0 is correct FASTA file."
	echo
	exit 1
	}
echo
echo "Preliminary probe sequences saved as $PROBEPRELIM for possible later usage."
confirmgo

# Step 9: Make the final quality control of the probe sequences

echo "Step 9 of the pipeline - removal of probe sequences sharing ≥90% sequence similarity"
echo

# Check for sequence similarity between the developed probe sequences with CD-HIT-EST

# Clustering exons with 100% sequence identity: retaining unclustered exons and, in case of 100% sequence identity, retaining the longest exon
echo "Checking sequence similarity between the probe sequences (exons) Detecting identical probe sequences and retaining the longest one in such a case. Retaining also the unclustered probe sequences."
echo
cd-hit-est -i $PROBEPRELIM -o $PROBEPRELIMCLUSTER100 -d 0 -c 1.0 -p 1 || {
	echo
	echo "Error! Checking of the sequence similarity failed. Aborting. Check if file $PROBEPRELIM is correct FASTA file."
	echo
	exit 1
	}
echo
echo "Clustered exons with 100% sequence identity were saved as $PROBEPRELIMCLUSTER100 for possible later usage."
confirmgo

# Clustering and removing exons with more than a certain sequence similarity
echo "Detecting and removing probe sequences (exons) that are similar to each other above a certain threshold"
echo
cd-hit-est -i $PROBEPRELIMCLUSTER100 -o $PROBEPRELIMCLUSTER90 -d 0 -c $CDHITSIM -p 1 -g 1 || {
	echo
	echo "Error! Checking of the probe sequence failed. Aborting. Check if file $PROBEPRELIMCLUSTER100 is correct FASTA file."
	echo
	exit 1
	}
echo
# Finding those clusters from a CD-HIT CLSTR file that include only one sequence or multiple sequences with 100% identity (in which case the longest sequence is choosen)
python $SCRIPTDIR/grab_singleton_clusters.py -i $PROBEPRELIMCLUSTER90.clstr -o $UNIQUEPROBEPRELIMCLUSTER90 || {
	echo
	echo "Error! Checking of the probe sequence failed. Aborting. Check if file $PROBEPRELIMCLUSTER90.clstr is correct output of cd-hit-est."
	echo
	exit 1
	}

echo "Unclustered exons and clustered exons with 100% identity were saved as $UNIQUEPROBEPRELIMCLUSTER90 for possible later usage."
confirmgo

echo "Postprocessing extracted sequences"
{ grep -v '>Cluster' $UNIQUEPROBEPRELIMCLUSTER90 | cut -d ' ' -f 2 | sed -e 's/\.\.\./\\\>/' -e 's/^/^/' > $UNIQUEPROBEPRELIM &&
grep -A 1 -f $UNIQUEPROBEPRELIM $PROBEPRELIMCLUSTER100 | sed '/^--$/d' > $UNIQUEPROBEPRELIMF; } || {
	echo
	echo "Error! Postprocessing of extracted sequences failed. Aborting. Check if file $UNIQUEPROBEPRELIMCLUSTER90 is correct output of cd-hit-est."
	echo
	exit 1
	}
echo
echo "Postprocessed extracted sequences were saved as $UNIQUEPROBEPRELIMF for possible later usage."
confirmgo

# Step 10: Retention of those probe sequences that comprise exons of a certain minimum length (default is 120 bp) and have a certain minimum total locus length

echo "Step 10 of the pipeline - retention of those probe sequences that comprise exons of a certain minimum length ($BAITL bp) and have a certain minimum total locus length"
echo

# One of the three outfiles of last steps of part 9 is a FASTA file, it has to be converted to TAB
echo "Converting FASTA from step 9 to TAB"
fasta2tab $UNIQUEPROBEPRELIMF $PROBEPRELIMCDHIT

echo "Writing the exons into temporal file"
{ awk '{print $1"\t"length($2)}' $PROBEPRELIMCDHIT | sed 's/_/\t/g' | cut -f 2,6 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' > $PROBEPRELIMCDHIT2; } || {
	echo
	echo "Error! Writing of the exons failed. Aborting. Check if file $PROBEPRELIMCDHIT is correct file (sequences in TAB)."
	echo
	exit 1
	}
echo

# Extract and sort the exons making up genes of a minimum total length
echo "Extracting and sorting the exons making up genes of ≥$MINLOCUSLENGTH bp"
{ sed 's/^/Assembly_/' $PROBEPRELIMCDHIT2 | cut -f 1 -d " " | sort -k 1,1 > $PROBEPRELIMFORJOIN; } || {
	echo
	echo "Error! Extraction and sorting of exons failed. Aborting. Check if file $PROBEPRELIMCDHIT2 is correct file (sequences in TAB)."
	echo
	exit 1
	}
echo

# Modify the exon number and sort
echo "Modifying the exon number and sorting"
{ sed 's/_C/\tC/' $PROBEPRELIMCDHIT | sort -k 1,1 > $PROBEPRELIMSORT; } || {
	echo
	echo "Error! Modification of the exons failed. Aborting. Check if file $PROBEPRELIMCDHIT is correct file (sequences in TAB)."
	echo
	exit 1
	}
echo

# Make a file with all exons of a certain minimum length making up genes of a certain minimum total length
echo "Joining all exons ≥$BAITL bp and making up genes of ≥$MINLOCUSLENGTH bp"
join $PROBEPRELIMFORJOIN $PROBEPRELIMSORT | sed 's/^\(.\+\) \(Contig\)/>\1_\2/' > $PROBEPRELIMFIN || {
	echo
	echo "Error! Joining of exons failed. Aborting. Check if files $PROBEPRELIMFORJOIN and $PROBEPRELIMSORT are correct files (sequences in TAB)."
	echo
	exit 1
	}
echo

# Probe design summary

# Calculation of the total number of base pairs
echo "Calculating the total number of base pairs"
echo "  $(awk '{print $1"\t"length($3)}' $PROBEPRELIMFIN | sed 's/_/\t/g' | cut -f 2,3 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' | awk '{s+=$3;c++}END{print s}') bp make up genes of ≥$MINLOCUSLENGTH bp."
confirmgo

# Calculation of the total number of genes
echo "Calculating the total number of genes"
echo "  There are $(awk '{print $1"\t"length($3)}' $PROBEPRELIMFIN | sed 's/_/\t/g' | cut -f 2,3 | awk '{a[$1]++;b[$1]+=$2}END{for (i in a) print i,a[i],b[i]}' | awk '$3>'"$MINLOCUSLENGTHN"'' | wc -l) genes in total."
confirmgo

# Calculation of the total number of exons
echo "Calculating the total number of exons"
echo "  Total number of exons ≥$BAITL bp: $(wc -l $PROBEPRELIMFIN | cut -f 1 -d " ")."
confirmgo

# Convert TAB to FASTA
echo "Converting TAB to FASTA"
{ sed 's/ /\n/' $PROBEPRELIMFIN > $PROBESEQUENCES; } || {
	echo
	echo "Error! Conversion of TAB to FASTA failed. Aborting. Check if file $PROBEPRELIMFIN is correct file (sequences in TAB)."
	echo
	exit 1
	}
echo

echo "Success!"
echo
echo "================================================================================"
echo "Probes with all sequences (including putative plastid genes) are in $PROBESEQUENCES"
echo "This file contains the probe sequences. In next step, putative plastid sequences will be removed. We STRONGLY RECOMMEND to remove those genes from the final probe set."
echo "================================================================================"
confirmgo

# Step 11 - removal of putative cpDNA sequences in final probe list

echo "Step 11 of the pipeline - detection of probe sequences sharing ≥$BLATIDENT% sequence imilarity with the plastome reference"
echo

# Remove remaining cp genes from probe set
echo "Detecting remaining plastid genes in probe set"
blat -t=dna -q=dna -minIdentity=$BLATIDENT -out=pslx $REFERENCECP $PROBESEQUENCES $PROBESEQUENCESCP || {
	echo
	echo "Error! Detection of remaining plastid genes failed. Aborting. Check if files $REFERENCECP and $PROBESEQUENCES are correct FASTA file."
	echo
	exit 1
	}
echo

echo "================================================================================"
echo "File $PROBESEQUENCESCP contains putative plastid genes found in $PROBESEQUENCES set. We STRONGLY RECOMMEND to remove those genes from the final probe set."
echo "================================================================================"
confirmgo

# Extracting names of putative plastid sequence
echo "Preparing to remove putative plastid genes from final probe set"
sed 1,5d $PROBESEQUENCESCP | cut -f 10 | sort -u | sed 's/^/\\</g' | sed 's/$/\\>/g' > $PROBESEQUENCESCPLIST || {
	echo
	echo "Error! Processing of file containing remaining putative plastid sequences failed. Aborting. Check if file $PROBESEQUENCESCP correct PSLX file (BLAT output)."
	echo
	exit 1
	}
# Removing putative plastid sequence from all probes, converting to FASTA
echo "Removing remaining putative plastid genes from final probe set and converting to FASTA"
grep -v -f $PROBESEQUENCESCPLIST $PROBEPRELIMFIN | sed 's/ /\n/' > $PROBESEQUENCESNOCP || {
	echo
	echo "Error! Removal of putative plastid genes from final probe set failed. Aborting. Check if file $PROBEPRELIMFIN is correct file (sequences in TAB)."
	echo
	exit 1
	}
echo

echo "Success!"
echo
echo "================================================================================"
echo "Final output file was written as $PROBESEQUENCESNOCP"
echo "This file contains the probe sequences. Putative plastid genes were removed."
echo "================================================================================"
confirmgo

# Remove temporal files
echo "Removing unneeded temporal files"
rm $REFERENCECP $SEQUENCES $SEQUENCESTAB $SEQUENCESTABASSE $SEQUENCESPROBESLOCUSLENGTH $SEQUENCESPROBESLOCUSLENGTHFORJOIN $SEQUENCESTABASSEBAITL $SEQUENCESTABASSEBAITLSORT $SEQUENCESPROBES120600FIN $SEQUENCESPROBES120600MODIF $SEQUENCESPROBES120600ASSEM $SEQUENCESPROBES120600CONTIG $PROBEPRELIM0 $PROBEPRELIMCLUSTER90 $UNIQUEPROBEPRELIM $PROBEPRELIMCLUSTER100 $UNIQUEPROBEPRELIMF $PROBEPRELIMCDHIT $PROBEPRELIMFORJOIN $PROBEPRELIMSORT $PROBEPRELIMFIN $PROBESEQUENCESCPLIST || {
	echo
	echo "Error! Removal of temporal files failed. Remove following files manually: \"$REFERENCECP\", \"$SEQUENCES\", \"$SEQUENCESTAB\", \"$SEQUENCESTABASSE\", \"$SEQUENCESPROBESLOCUSLENGTH\",\"$SEQUENCESPROBESLOCUSLENGTHFORJOIN\", \"$SEQUENCESTABASSEBAITL\", \"$SEQUENCESTABASSEBAITLSORT\",\"$SEQUENCESPROBES120600FIN\", \"$SEQUENCESPROBES120600MODIF\", \"$SEQUENCESPROBES120600ASSEM\", \"$SEQUENCESPROBES120600CONTIG\", \"$PROBEPRELIM0\", \"$PROBEPRELIMCLUSTER90\", \"$UNIQUEPROBEPRELIM\", \"$PROBEPRELIMCDHIT\", \"$PROBEPRELIMFORJOIN\", \"$PROBEPRELIMSORT\", \"$PROBEPRELIMFIN\" and \"$PROBESEQUENCESCPLIST\"."
	confirmgo
	}
echo

# List kept files which user can use for another analysis
echo "Following files are kept for possible later usage (see PDF manual for details):"
echo "================================================================================"
echo "1) Preliminary probe sequences:"
echo "$PROBEPRELIM"
echo "2) Unclustered exons and clustered exons with 100% sequence identity:"
echo "$PROBEPRELIMCLUSTER100"
echo "3) Clustered exons with 100% sequence identity:"
echo "$UNIQUEPROBEPRELIMCLUSTER90"
echo "4) Unclustered exons / exons with less than a certain sequence similarity:"
echo "$UNIQUEPROBEPRELIMF"
echo "5) Contigs comprising exons ≥ bait length having a certain minimum total locus length:"
echo "$PROBEPRELIMCDHIT2"
echo "6) Probe sequences in FASTA format (including putative plastid genes):"
echo "$PROBESEQUENCES"
echo "7) Putative plastid genes in final probe set:"
echo "$PROBESEQUENCESCP"
echo "8) Final probe sequences in FASTA format:"
echo "$PROBESEQUENCESNOCP"
echo "================================================================================"
confirmgo

echo "Script exited successfully..."
echo

exit
