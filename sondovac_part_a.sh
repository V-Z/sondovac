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

echo "This is part A of the pipeline."
echo
echo "This part is for filtering of raw data and their preparation for assembly in"
echo "Geneious. Results of Geneious assembly are processed in part B to get the final"
echo "list of low-copy nuclear probe sequences. See README and/or manual for details."

# Default values
# Counter if not both -i and -n options are used
CHECKMODE=0
# If not specifying explicitly otherwise (using -n), running in interactive mode
STARTINI="I"
# flash -M maximum overlap length expected in approximately 90% of read pairs
FLASHM=65
# BLAT -minIdentity between the unique transcripts and the genome skim data
BLATIDENT=85
# Remove transcripts with >1000 BLAT scores (or another value selected by user)
BLATSCORE=1000
# Default name of output files
OUTPUTFILENAME=`realpath "output"`
# Usage of genome skim (22) or transcript sequences (23, parameter "-g")
PSLXCUT=22

# Create empty variables for file names
INPUTFILE=""
INPUTFILE0=""
INPUTFILE1=""
REFERENCECP=""
REFERENCECP0=""
INPUTFQ1=""
INPUTFQ2=""
REFERENCEMT=""
REFERENCEMT0=""

# Parse initial arguments
while getopts "hvulrpeo:inf:c:m:t:q:a:y:s:g" START; do
	case "$START" in
		h|v)
			generaloptions
			echo
			echo -e "\tIf options -f, -c, -m, -t and/or -q are used and the script is running"
			echo -e "\t  in interactive mode, those values will be used as defaults, but may"
			echo -e "\t  later be overwritten."
			echo
			echo -e "\tOptions required for running in non-interactive mode:"
			echo -e "\t-f\tTranscriptome input file in FASTA format."
			echo -e "\t-c\tPlastome reference sequence input file in FASTA format."
			echo -e "\t-m\tMitochondriome reference sequence input file in FASTA format."
			echo -e "\t\t  This file is optional. In interactive mode you will each time"
			echo -e "\t\t  be asked if you wish to use it."
			echo -e "\t-t\tPaired-end genome skim input file in FASTQ format (first file)."
			echo -e "\t-q\tPaired-end genome skim input file in FASTQ format (second file)."
			echo
			echo -e "\tOther optional arguments (if not provided, default values are used):"
			echo -e "\t-a\tMaximum overlap length expected in approximately 90% of read"
			echo -e "\t\t  pairs (parameter \"-M\" of FLASH, see its manual for details)."
			echo -e "\t\tDefault value: 65 (integer ranging from 10 to 300)"
			echo -e "\t-y\tSequence similarity between unique transcripts and the filtered,"
			echo -e "\t\t  combined genome skim reads (parameter \"-minIdentity\" of BLAT,"
			echo -e "\t\t  see its manual for details)."
			echo -e "\t\tDefault value: 85 (integer ranging from 70 to 100; the default"
			echo -e "\t\t  value of 85% minimum sequence similarity suggests gene"
			echo -e "\t\t  orthology)"
			echo -e "\t-s\tNumber of BLAT hits per transcript when matching unique"
			echo -e "\t\t  transcripts and the filtered, combined genome skim reads."
			echo -e "\t\tDefault value: 1000 (integer ranging from 100 to 10000)"
			echo -e "\t-g\tUse genome skim sequences instead of transcripts for making the"
			echo -e "\t\t  probes. Default is usage of genome skim sequences (no parameter)."
			echo -e "\tWARNING! If parameters -a, -y, -s or -g are not provided, default values"
			echo -e "\t\t are taken, and it is not possible to change them later (not even"
			echo -e "\t\t in interactive mode)."
			echo
			echo "Examples:"
			echo "Basic and the most simple usage:"
			echo "$0 -i"
			echo "Specify some of required input files, otherwise run interactively:"
			echo "$0 -i -f input.fa -t reads1.fastq -q reads2.fastq"
			echo "Running in non-interactive, automated mode:"
			echo "$0 -n -f input.fa -c referencecp.fa -m referencemt.fa -t reads1.fastq -q reads2.fastq"
			echo "Modify parameter -a, otherwise run interactively:"
			echo "$0 -i -a 300"
			echo "Run in non-interactive mode (parameter -n) - in such case user must specify all"
			echo "  required input files (parameters -f, -c, -m, -t and -q). Moreover, parameter"
			echo "  -y is modified:"
			echo "$0 -n -f input.fa -c referencecp.fa -m referencemt.fa -t reads1.fastq -q reads2.fastq -y 90"
			echo "Modifying parameter -s. Note interactive mode -i is implicit and does not need"
			echo "  to be specified explicitly:"
			echo "$0 -s 950"
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
		f)
			INPUTFILE1=$OPTARG
			echo "Transcriptome file: $INPUTFILE1"
			;;
		c)
			REFERENCECP0=$OPTARG
			echo "Plastome reference: $REFERENCECP0"
			;;
		m)
			REFERENCEMT0=$OPTARG
			echo "Mitochondriome reference: $REFERENCEMT0"
			;;
		t)
			INPUTFQ1=$OPTARG
			echo "FASTQ reads 1: $INPUTFQ1"
			;;
		q)
			INPUTFQ2=$OPTARG
			echo "FASTQ reads 2: $INPUTFQ2"
			;;
		a)
			FLASHM=$OPTARG
			# Check if provided value makes sense
			if [[ "$FLASHM" =~ ^[0-9]+$ ]] && [ "$FLASHM" -ge 10 -a "$FLASHM" -le 300 ]; then
				echo "Maximum overlap length expected in approximately 90% of read pairs: $FLASHM"
				else
					echo "Error! For parameter \"-a\" you did not provide an integer ranging from 10 to 300!"
					echo
					exit 1
					fi
			;;
		y)
			BLATIDENT=$OPTARG
			# Check if provided value makes sense
			if [[ "$BLATIDENT" =~ ^[0-9]+$ ]] && [ "$BLATIDENT" -ge 70 -a "$BLATIDENT" -le 100 ]; then
				echo "BLAT score for identity between unique transcripts and genome skim data: $BLATIDENT"
				else
					echo
					echo "Error! For parameter \"-y\" you did not provide an integer of range from 70 to 100!"
					echo
					exit 1
					fi
			;;
		s)
			BLATSCORE=$OPTARG
			# Check if provided value makes sense
			if [[ "$BLATSCORE" =~ ^[0-9]+$ ]] && [ "$BLATSCORE" -ge 100 -a "$BLATSCORE" -le 10000 ]; then
				echo "BLAT score: $BLATSCORE"
				else
					echo
					echo "Error! For parameter \"-s\" you did not provide an integer ranging from 100 to 10000!"
					echo
					exit 1
					fi
			;;
		g)
			PSLXCUT=23
			echo "The script will use transcript sequences instead of genome skim sequences."
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

# Ensure user reads introductory information
confirmgo

# Check operating system
oscheck

# Set variables for working directory and PATH
workdirpath

# Check availability of all needed binaries

# Check if realpath is available
checktools realpath

# Check if paste is available
checktools paste

# Check if cut is available
checktools cut

# Check if uniq is available
checktools uniq

# Check if awk is available
checktools awk

# Check if sort is available
checktools sort

# Check if join is available
checktools join

# Check if sed is available
checktools sed

# Check if grep is available
checktools grep

# Check if cat is available
checktools cat

# Check if perl is available
checktools perl

# Check if ls is available
checktools ls

# Check if cp is available
checktools cp

# Check if mkdir is available
checktools mkdir

# Check if tr is available
checktools tr

# Check if BLAT is available
checkblat

# Function to compile Bowtie2
BOWTIE2V='2.3.3.1'
function compilebowtie {
	{
	echo &&
	checktools make &&
	checktools g++ &&
	echo &&
	echo "Compiling Bowtie2 from source code..." &&
	echo &&
	cd $1 &&
	make -s &&
	cp -r bowtie2* $BIN/ &&
	cd $WORKDIR &&
	echo &&
	echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\"" &&
	echo "   and \"bowtie2-build-s\" are available. OK."
	} || {
		echo
		echo "Error! Compilation failed. Please, go to"
		echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/"
		echo "  download bowtie2-$BOWTIE2V-source.zip, compile it, and ensure it is in PATH."
		echo
		exit 1
		}
	}

# Check if Bowtie2 is available
{ command -v bowtie2 >/dev/null 2>&1 &&
	command -v bowtie2-align-l >/dev/null 2>&1 &&
	command -v bowtie2-align-s >/dev/null 2>&1 &&
	command -v bowtie2-build >/dev/null 2>&1 &&
	command -v bowtie2-build-l >/dev/null 2>&1 &&
	command -v bowtie2-build-s >/dev/null 2>&1 &&
	echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\"" &&
	echo "  and \"bowtie2-build-s\" are available. OK."
	} || {
		echo
		echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\""
		echo "  and \"bowtie2-build-s\" are required but not installed or available in PATH."
		if [ "$STARTINI" == "I" ]; then
			echo
			echo "Type \"C\" to compile \"Bowtie2-$BOWTIE2V\" from source available together with this script."
			echo "Type \"S\" to compile \"Bowtie2-$BOWTIE2V\" from source code downloaded from"
			echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/"
			echo "Type \"D\" to download \"Bowtie2-$BOWTIE2V\" binary from"
			echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/ for your OS."
			echo "Type \"B\" to copy Bowtie2-2.$BOWTIE2V binary available together with the script"
			echo "  (recommended, available for Linux and Mac OS X)."
			echo "Type \"H\" for installation using Homebrew (only for Mac OS X, recommended)."
			echo "  See \"brew info homebrew/science/bowtie2\" for more details."
			echo "Type \"M\" for manual installation - script will exit, and you will have to install"
			echo "  Bowtie2 yourself."
			read BOWTIE
			while :
				do
					case "$BOWTIE" in
						C|c)
							compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
							break
							;;
						D|d)
							if [ "$OS" == "Linux" ]; then
								{
								echo "Downloading \"Bowtie2-$BOWTIE2V\" binaries for $OS $OSB" &&
								$DOWNLOADER bowtie2-$BOWTIE2V-linux-x86_64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-linux-x86_64.zip &&
								unzip -nq bowtie2-$BOWTIE2V-linux-x86_64.zip &&
								cp bowtie2-$BOWTIE2V/bowtie2* $BIN/ &&
								echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\"" &&
								echo "  and \"bowtie2-build-s\" are available. OK."
								} || {
									echo
									echo "Error! Download failed. Please, go to"
									echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/"
									echo "  and download and unpack bowtie2-$BOWTIE2V-linux-x86_64.zip manually."
									echo
									exit 1
									}
								elif [ "$OS" == "Mac" ]; then
									{
									echo "Downloading \"Bowtie2-$BOWTIE2V\" binaries for $OS $OSB" &&
									$DOWNLOADER bowtie2-$BOWTIE2V-macos-x86_64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-macos-x86_64.zip &&
									unzip -nq bowtie2-$BOWTIE2V-macos-x86_64.zip &&
									cp bowtie2-$BOWTIE2V/bowtie2* $BIN/ &&
									echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\"" &&
									echo "  and \"bowtie2-build-s\" are available. OK."
									} || {
										echo
										echo "Error! Download failed. Please, go to"
										echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/"
										echo "  and download and unpack bowtie2-$BOWTIE2V-macos-x86_64.zip manually."
										echo
										exit 1
										}
									elif [ "$OS" == "Windows" ]; then
										{
										echo "Downloading \"Bowtie2-$BOWTIE2V\" binaries for $OS $OSB" &&
										$DOWNLOADER bowtie2-$BOWTIE2V-mingw-win64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-mingw-win64.zip &&
										unzip -nq bowtie2-$BOWTIE2V-mingw-win64.zip &&
										cp bowtie2-$BOWTIE2V/bowtie2* $BIN/ &&
										echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\""
										echo "  and \"bowtie2-build-s\" are available. OK."
										} || {
											echo
											echo "Error! Download failed. Please, go to"
											echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/"
											echo "  and download and unpack bowtie2-$BOWTIE2V-mingw-win64.zip manually."
											echo
											exit 1
											}
										else
											echo "Unknown OS or OS without \"Bowtie2-$BOWTIE2V\" binary available."
											compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
											fi
							break
							;;
						S|s)
							echo
							downloaderselector
							checktools unzip
							echo
							echo "Downloading \"Bowtie2-$BOWTIE2V\" source code"
							$DOWNLOADER bowtie2-$BOWTIE2V-source.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-source.zip || {
								echo
								echo "Error! Download failed. Please, go to"
								echo "  http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/"
								echo "  and download bowtie2-$BOWTIE2V-source.zip and compile it manually."
								echo
								exit 1
								}
							unzip -nq $WORKDIR/bowtie2-$BOWTIE2V-source.zip
							compilebowtie bowtie2-$BOWTIE2V
							break
							;;
						B|b)
							case "$OS" in
								Mac)
									echo "Copying Bowtie-$BOWTIE2V binaries"
									cp -pr $SCRIPTDIR/pkgs/macosx/bin/bowtie2* $BIN/
									mkdir -p $WORKDIR/bin/share
									cp -pr $SCRIPTDIR/pkgs/macosx/share/man $WORKDIR/bin/share/
									break
									;;
								Linux)
									echo "Copying Bowtie-$BOWTIE2V binaries"
									cp -pr $SCRIPTDIR/pkgs/linux64b/bin/bowtie2* $BIN/
									mkdir -p $WORKDIR/bin/share
									cp -pr $SCRIPTDIR/pkgs/linux64b/share/man $WORKDIR/bin/share/
									break
									;;
								*) echo
									echo "Binary is not available for $OS $OSB. Going to compile it from source code."
									echo
									compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
									;;
								esac
							break
							;;
						H|h)
							if [ "$OS" == "Mac" ]; then
								{ echo "Installing \"Bowtie2-$BOWTIE2V\" using Homebrew" &&
								brew install homebrew/science/bowtie2 &&
								echo "\"Bowtie2\" is available. OK."
								} || {
									echo
									echo "Error! Installation of \"Bowtie2\" failed. Please, do it manually. For details see"
									echo "\"brew info homebrew/science/bowtie2\" and \"brew help\"."
									echo
									exit 1
									}
								else
									echo "This is not Mac OS X. Going to compile..."
									compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
									fi
							break
							;;
						M|m)
							echo "Please, go to http://bowtie-bio.sourceforge.net/bowtie2/index.shtml and install"
							echo " latest Bowtie2, and ensure it is in PATH."
							exit 2
							;;
						*) echo "Wrong option. Use C, D, S, B, H or M." && read BOWTIE;;
						esac
					done
					else
						exit 1
						fi
		}

function compilesamtools {
	cd $SCRIPTDIR/src/samtools-1.6 &&
	echo &&
	checktools make &&
	echo &&
	echo "Compiling samtools..." &&
	./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.6 --prefix=$WORKDIR/bin --without-curses &&
	make -s all plugins-htslib all-htslib &&
	make -s install install-htslib &&
	cd $WORKDIR &&
	echo "\"samtools\" is available. OK."
	}

# Check if SAMtools is available
{ command -v samtools >/dev/null 2>&1 && echo "\"samtools\" is available. OK."; } || {
	echo
	echo >&2 "\"samtools\" is required but not installed or available in PATH."
	if [ "$STARTINI" == "I" ]; then
		echo
		echo "Type \"C\" to compile \"SAMtools-1.6\" from source available together with this script."
		echo "  Makefile was modified not to require GNU ncurses library."
		echo "Type \"S\" to download latest developmental \"SAMtools\" source from"
		echo "  https://github.com/samtools/samtools/ and compile it. Compilation requires GNU"
		echo "  ncurses library and is recommended only for advanced users. If compilation"
		echo "  fails, check SAMtools' INSTALL file for details and adjust its Makefile."
		echo "Type \"B\" to copy \"SAMtools-1.6\" binary available together with the script"
		echo "  (recommended, available for Linux and Mac OS X)."
		echo "Type \"H\" for installation using Homebrew (only for Mac OS X, recommended)."
		echo "  See \"brew info homebrew/science/samtools\" for more details."
		echo "Type \"M\" for manual installation - script will exit, and you will have to"
		echo "  install \"SAMtools\" yourself."
		read SAMTOOLS
		while :
			do
				case "$SAMTOOLS" in
					C|c)
						compilesamtools || {
							echo
							echo "Error! Compilation failed. Please, go to http://www.htslib.org/download/"
							echo "  download samtools-1.6, compile it, and ensure it is in PATH."
							echo
							exit 1
							}
						break
						;;
					S|s)
						{
						echo
						downloaderselector &&
						checktools unzip &&
						checktools pwd &&
						checktools make &&
						checktools gcc &&
						echo &&
						echo "Downloading \"SAMtools\" sources..." &&
						$DOWNLOADER develop.zip https://github.com/samtools/samtools/archive/develop.zip &&
						unzip -nq develop.zip &&
						cd samtools-develop &&
						$DOWNLOADER develop.zip https://github.com/samtools/htslib/archive/develop.zip &&
						unzip -nq develop.zip &&
						echo "Compiling SAMtools" &&
						make -s HTSDIR=`pwd`/htslib-develop &&
						make -s prefix=$WORKDIR/bin install &&
						cd $WORKDIR &&
						echo "\"samtools\" is available. OK."
						} || {
							echo
							echo "Error! Compilation failed. Please, go to http://www.htslib.org/download/ for"
							echo "  latest stable version of \"SAMtools\" or https://github.com/samtools/ for latest"
							echo "  developmental version, download samtools, compile it, and ensure it is in PATH."
							echo
							exit 1
							}
						break
						;;
					B|b)
						case "$OS" in
							Mac)
								echo "Copying \"SAMtools\" binaries"
								cp -p $SCRIPTDIR/pkgs/macosx/bin/ace2sam $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/blast2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/bowtie2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/export2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/interpolate_sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/maq2sam-long $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/maq2sam-short $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/md5fa $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/md5sum-lite $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/novo2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/plot-bamstats $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/psl2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/samtools $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/samtools.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/sam2vcf.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/seq_cache_populate.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/soap2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/varfilter.py $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/wgsim $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/wgsim_eval.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/macosx/bin/zoom2sam.pl $BIN/
								mkdir -p $WORKDIR/bin/share/man/man1
								cp -p $SCRIPTDIR/pkgs/macosx/share/man/man1/samtools.1 $WORKDIR/bin/share/man/man1/
								break
								;;
							Linux)
								echo "Copying \"SAMtools\" binaries"
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/ace2sam $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/blast2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/bowtie2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/export2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/interpolate_sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/maq2sam-long $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/maq2sam-short $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/md5fa $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/md5sum-lite $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/novo2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/plot-bamstats $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/psl2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/samtools $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/samtools.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/sam2vcf.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/seq_cache_populate.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/soap2sam.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/varfilter.py $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/wgsim $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/wgsim_eval.pl $BIN/
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/zoom2sam.pl $BIN/
								mkdir -p $WORKDIR/bin/share/man/man1
								cp -p $SCRIPTDIR/pkgs/linux64b/share/man/man1/samtools.1 $WORKDIR/bin/share/man/man1/
								;;
							*) echo
								echo "Binary is not available for $OS $OSB."
								echo
								compilesamtools
								;;
							esac
						break
						;;
					H|h)
						if [ "$OS" == "Mac" ]; then
						{ echo "Installing SAMtools using Homebrew" &&
						brew install homebrew/science/samtools &&
						echo "\"SAMtools\" is available. OK."
						} || {
							echo
							echo "Error! Installation of \"SAMtools\" failed. Please, do it manually. For details see"
							echo "\"brew info homebrew/science/samtools\" and \"brew help\"."
							echo
							exit 1
							}
							else
								echo "This is not Mac OS X. Going to compile..."
								compilesamtools
								fi
						break
						;;
					M|m)
						echo "Please, go to http://www.htslib.org/ and install \"SAMtools\" and ensure it is"
						echo "  in PATH."
						exit 2
						;;
					*) echo "Wrong option. Use C, S, B, H or M." && read SAMTOOLS;;
					esac
				done
			else
				exit 1
				fi
	}

# Function to compile FLASH
function compileflash {
	{
	echo &&
	checktools make &&
	checktools gcc &&
	echo &&
	echo "Compiling \"FLASH\" from source code..." &&
	echo &&
	cd $1 &&
	make -s &&
	cp flash $BIN/ &&
	cd $WORKDIR &&
	echo &&
	echo "\"flash\" is available. OK."
	} || {
		echo
		echo "Error! Compilation failed. Please, go to"
		echo "  http://sourceforge.net/projects/flashpage/files/ download latest"
		echo "  FLASH-*.tar.gz, compile it manually, and ensure it is in PATH."
		echo
		exit 1
		}
	}

# Check if FLASH is available
{ command -v flash >/dev/null 2>&1 && echo "\"flash\" is available. OK."; } || {
	echo
	echo >&2 "FLASH is required but not installed or available in PATH."
	if [ "$STARTINI" == "I" ]; then
		echo
		echo "Type \"C\" to compile \"FLASH\" from source available together with this script."
		echo "Type \"S\" to download \"FLASH\" source from"
		echo "  http://sourceforge.net/projects/flashpage/ and compile it."
		echo "Type \"D\" to download \"FLASH\" binary from"
		echo "  http://sourceforge.net/projects/flashpage/ (available only for Windows)."
		echo "Type \"B\" to copy \"FLASH 1.2.11\" binary available together with the script"
		echo "  (recommended, available for Linux and Mac OS X)."
		echo "Type \"H\" for installation using Homebrew (only for Mac OS X, recommended)."
		echo "  See \"brew info homebrew/science/flash\" for more details."
		echo "Type \"M\" for manual installation - script will exit, and you will have to"
		echo "  install \"FLASH\" yourself."
		read FLASH
		while :
			do
				case "$FLASH" in
					C|c)
						compileflash $SCRIPTDIR/src/FLASH-1.2.11
						break
						;;
					S|s)
						echo
						downloaderselector
						checktools tar
						checktools gunzip
						echo
						echo "Downloading \"FLASH\" source code"
						$DOWNLOADER FLASH-1.2.11.tar.gz http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz || {
							echo
							echo "Error! Download failed. Please, go to http://sourceforge.net/projects/flashpage/files/"
							echo "  download FLASH-*.tar.gz, compile it manually, and ensure it is in PATH."
							echo
							exit 1
							}
						tar xzf FLASH-1.2.11.tar.gz
						compileflash FLASH-1.2.11
						break
						;;
					D|d)
						if [ "$OS" != "Windows" ]; then
							echo "OS is not Windows!"
							compileflash $SCRIPTDIR/src/FLASH-1.2.11
							else
								echo
								downloaderselector
								checktools unzip
								checktools chmod
								checktools mv
								echo
								echo "Downloading \"FLASH\" for $OS"
								$DOWNLOADER FLASH-1.2.11-windows-bin.zip http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11-windows-bin.zip || {
									echo
									echo "Error! Download failed. Please, go to http://sourceforge.net/projects/flashpage/files/"
									echo "  download FLASH-*windows-bin.zip, unpack it, and ensure it is in PATH."
									echo
									exit 1
									}
								unzip -nq FLASH-1.2.11-windows-bin.zip
								chmod +x flash.exe
								mv flash.exe $BIN/
								fi
						break
						;;
					B|b)
						case "$OS" in
							Mac)
								echo "Copying \"FLASH\" binary"
								cp -p $SCRIPTDIR/pkgs/macosx/bin/flash $BIN/
								break
								;;
							Linux)
								echo "Copying \"FLASH\" binary"
								cp -p $SCRIPTDIR/pkgs/linux64b/bin/flash $BIN/
								break
								;;
							*) echo
								echo "Binary is not available for $OS $OSB."
								echo
								compileflash $SCRIPTDIR/src/FLASH-1.2.11
								;;
							esac
						break
						;;
					H|h)
						if [ "$OS" == "Mac" ]; then			
							{ echo "Installing \"FLASH\" using Homebrew" &&
							brew install homebrew/science/flash &&
							echo "\"FLASH\" is available. OK."
							} || {
								echo
								echo "Error! Installation of \"FLASH\" failed. Please, do it manually. For details see"
								echo "\"brew info homebrew/science/flash\" and \"brew help\"."
								echo
								exit 1
								}
							else
								echo "This is not Mac OS X. Going to compile..."
								compileflash $SCRIPTDIR/src/FLASH-1.2.11
								fi
						break
						;;
					M|m)
						echo "Please, go to http://ccb.jhu.edu/software/FLASH/ and install \"FLASH\" and ensure it"
						echo "  is in PATH"
						exit 2
						;;
					*) echo "Wrong option. Use C, S, D, B, H or M." && read FLASH;;
					esac
				done
		else
			exit 1
			fi
	}

echo
# Input files
CHECKFILEREADOUT=""

# Input data, transcriptome data in FASTA format
readinputfile -f "transcriptome input file in FASTA format" $INPUTFILE1
INPUTFILE1=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, plastome reference in FASTA format
readinputfile -c "plastome reference sequence input file in FASTA format" $REFERENCECP0
REFERENCECP0=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, genome skim reads in FASTQ, file 1
readinputfile -m "paired-end genome skim input file in FASTQ format (first file)" $INPUTFQ1
INPUTFQ1=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, genome skim reads in FASTQ, file 2
readinputfile -t "paired-end genome skim input file in FASTQ format (second file)" $INPUTFQ2
INPUTFQ2=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, mitochondriome reference in FASTA
if [ "$STARTINI" == "I" ]; then
	echo
	echo "Would you like to use mitochondriome reference sequence input file in FASTA"
	echo "  format? (Yes/No)"
	read MTINPUTQ
	while :
	do
		case $MTINPUTQ in
			Y|y|Yes|yes|YES)
				readinputfile -q "mitochondriome reference sequence input file in FASTA format" $REFERENCEMT0
				REFERENCEMT0=$CHECKFILEREADOUT
				CHECKFILEREADOUT=""
				break
				;;
			N|n|No|no|NO)
				echo
				echo "OK, we will not use mitochondriome reference sequence. Continuing."
				REFERENCEMT0=""
				REFERENCEMT=""
				break
				;;
			*) echo "Wrong option. Use Y or N." && read $MTINPUTQ;;
			esac
		done
	echo
	fi

# Notify user if mitochondriome is missing
if [ -z "$REFERENCEMT0" ]; then
	echo
	echo "Warning! There is no mitochondriome reference sequence."
	echo
	REFERENCEMT=""
	REFERENCEMT0=""
	else
		# Test if input file is readable
		if [[ -f $REFERENCEMT0 && -r $REFERENCEMT0 && -s $REFERENCEMT0 ]]; then
			echo
			echo "Input file \"$REFERENCEMT0\" exists and is readable. Proceeding..."
			else
				echo "Error! File \"$REFERENCEMT0\" does not exist, is empty or is not readable!"
				echo
				exit 1
				fi
		echo
		fi

# Input file in FASTA format
echo "Input file: $INPUTFILE1"
# Input file in FASTA format - checked not to be interleaved - temporary file - will be deleted
INPUTFILE0="${INPUTFILE1%.*}_non-interleaved.fasta"
# Input file in FASTA format - checked and renamed labels
INPUTFILE="${INPUTFILE1%.*}_renamed.fasta"
# List of old and new names of the transcriptome FASTA sequences
TRANSCRIPTOMEFASTANAMES="${INPUTFILE1%.*}_old_and_new_names.tsv"
# Output of BLAT (removal of transcripts sharing ≥90% sequence similarity)
BLATOUT="${OUTPUTFILENAME%.*}_blat_unique_transcripts.psl"
# List of unique transcripts - temporary file - will be deleted
UNIQUELIST="${OUTPUTFILENAME%.*}_trans-trans_unique_transcripts_sorted.txt"
# Input file converted into TXT - temporary file - will be deleted
INPUTTAB="${OUTPUTFILENAME%.*}.txt"
# Sorted input file in TXT - temporary file - will be deleted
SORTEDINPUT="${OUTPUTFILENAME%.*}_sorted.txt"
# Joined unique transcripts - temporary file - will be deleted
JOINEDTS="${OUTPUTFILENAME%.*}_unique_transcripts_trans-trans_plus_sequence.txt"
# Joined unique transcripts in tabular format - temporary file - will be deleted
JOINEDTABS="${OUTPUTFILENAME%.*}_tabs.txt"
# Joined unique transcripts in FASTA format
JOINEDFA="${OUTPUTFILENAME%.*}_unique_transcripts.fasta"
# Input - reference genome - cpDNA
echo "Input file: $REFERENCECP0"
# Input - reference genome - cpDNA - temporary file - will be deleted
REFERENCECP="${REFERENCECP0%.*}_non-interleaved.fasta"
# Reference genome - plastome index - temporary file - will be deleted
REFERENCECP2="${OUTPUTFILENAME%.*}.cp"
# Input reads in FASTQ
echo "Input file: $INPUTFQ1"
echo "Input file: $INPUTFQ2"
# cpDNA reads mapped to reference - temporary file - will be deleted
BOWTIE2CP="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads.sam"
# Genome skim data without cpDNA reads
FASTQNOCP="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads"
# Input - reference genome - mtDNA
if [ -n "$REFERENCEMT0" ]; then
	echo "Input file: $REFERENCEMT0"
	# Input - reference genome - mtDNA - temporary file - will be deleted
	REFERENCEMT="${REFERENCEMT0%.*}_non-interleaved.fasta"
	# Reference genome - mitochondriome index - temporary file - will be deleted
	REFERENCEMT2="${OUTPUTFILENAME%.*}.mt"
	# mtDNA reads mapped to reference - temporary file - will be deleted
	BOWTIE2MT="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_no_mt_reads.sam"
	# Genome skim data without mtDNA reads
	FASTQNOMT="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_no_mt_reads"
	# Combined paired-end genome skim reads
	fi
FLASHOUT="${OUTPUTFILENAME%.*}_combined_reads_no_cp_no_mt_reads"
# Output of BLAT (matching of the unique transcripts and the filtered, combined genome skim reads sharing ≥85% sequence similarity)
BLATOUTFIN="${OUTPUTFILENAME%.*}_blat_unique_transcripts_versus_genome_skim_data.pslx"
# Matching sequences in FASTA
BLATOUTFIN2="${OUTPUTFILENAME%.*}_blat_unique_transcripts_versus_genome_skim_data.fasta"
# FASTA converted into TSV - temporary file - will be deleted
TAB="${OUTPUTFILENAME%.*}_final.tab"
# Number of times each transcript hit a genome skim read - will be deleted
TABLIST="${OUTPUTFILENAME%.*}_transcript_hits.txt"
# Listed transcripts with >1000 BLAT score - will be deleted
TABBLAT="${OUTPUTFILENAME%.*}_1k_transcripts"
# Transcripts without >1000 BLAT score - will be deleted
TABREMOVED="${OUTPUTFILENAME%.*}_1k_transcripts-removed.tab"
# Final FASTA sequences for usage in Geneious
FINALA="${OUTPUTFILENAME%.*}_blat_unique_transcripts_versus_genome_skim_data-no_missing_fin.fsa"

# Check EOL of input files
echo
eolcheck $INPUTFILE1
eolcheck $REFERENCECP0
eolcheck $INPUTFQ1
eolcheck $INPUTFQ2
if [ -n "$REFERENCEMT0" ]; then
	eolcheck $REFERENCEMT0
	fi

# Check if FASTA input files are non-interleaved (required) - if not, FASTA input file converted
echo
echo "Checking if input FASTA files are non-interleaved (required) - interleaved"
echo "  FASTA files are converted not to be interleaved"
echo
noninterleavedfasta $INPUTFILE1 $INPUTFILE0
noninterleavedfasta $REFERENCECP0 $REFERENCECP
if [ -n "$REFERENCEMT0" ]; then
	noninterleavedfasta $REFERENCEMT0 $REFERENCEMT
	fi

# Transcriptome input file has required labeling scheme - only unique numbers

echo "FASTA sequence names in input file $INPUTFILE0"
echo "  must be renamed to be correctly handled in part A. New file with correct labels"
echo "  (only increasing unique numbers) will be created."
echo "Depending on size of your transcriptome file it can take longer time."
while read LINE; do
	N=$((++N)) &&
	echo $LINE | sed -e 's/^>.*$/>'$N'/'
	done < $INPUTFILE0 > $INPUTFILE || {
		echo "Error! Renaming failed. Aborting. Ensure $INPUTFILE0"
		echo "  has as labels of sequences only unique numbers (nothing else)."
		echo
		exit 1
		}
# Giving to user list of old and new labels
grep '^>' $INPUTFILE0 > transcript_fasta_labels_old
grep '^>' $INPUTFILE > transcript_fasta_labels_new
echo -e "Old FASTA labels\tNew FASTA labels" > $TRANSCRIPTOMEFASTANAMES
paste transcript_fasta_labels_old transcript_fasta_labels_new >> $TRANSCRIPTOMEFASTANAMES
sed -i 's/>//g' $TRANSCRIPTOMEFASTANAMES
rm transcript_fasta_labels_old transcript_fasta_labels_new
echo
echo "File $TRANSCRIPTOMEFASTANAMES contains two columns:"
echo "  1) Old labels in original $INPUTFILE0 and"
echo "  2) New labels in required format as in $INPUTFILE"
echo "  The sequences remain intact. This might be needed to trace back some sequences."
confirmgo

# Step 1: Obtain unique transcripts.

echo
echo "Step 1 of the pipeline - removal of transcripts sharing ≥90% sequence similarity."

# BLAT between the transcriptome itself
echo
echo "Launching BLAT between the transcriptome itself"
blat -t=dna -q=dna -noHead -out=psl $INPUTFILE $INPUTFILE $BLATOUT || {
	echo
	echo "Error! BLAT failed. Aborting. Check if $INPUTFILE is correct."
	echo
	exit 1
	}

# Filtering for unique transcripts
echo
echo "Filtering for unique transcripts:"
{ cut -f 10 $BLATOUT | uniq -c | awk '{print$1}' | sort | uniq -c; } || {
	echo
	echo "Error! Filtering of BLAT output failed. Aborting."
	echo "Check if files $INPUTFILE and $BLATOUT are correct."
	echo
	exit 1
	}
echo
echo "Filtered transcripts saved  for possible later usage as"
echo "  $BLATOUT  for possible later usage."
confirmgo

# Make a list of these unique transcripts (names and sequences) and convert this file to FASTA
echo "Making list of unique transcripts"
cut -f 10 $BLATOUT | uniq -c | awk '{if($1==1){print $0}}' | awk '{print $2}' | awk '{printf "%012d\n", $0;}' > $UNIQUELIST || {
	echo
	echo "Error! Making list of unique transcripts failed. Aborting."
	echo "Check if files $BLATOUT and $INPUTFILE are correct."
	echo
	exit 1
	}

# In order to use the join command, the original transcriptome file has to be converted to TXT, the transcript numbers have to be adjusted and the file sorted
echo
echo "Converting original data into TXT for subsequent joining"
fasta2tab $INPUTFILE $INPUTTAB || {
	echo
	echo "Error! Conversion of $INPUTFILE failed. Aborting."
	echo "Check if $INPUTFILE is valid FASTA file."
	echo
	exit 1
	}
echo

echo "Sorting unique transcripts"
{ awk '{$1=sprintf("%012d", $1); print $0}' $INPUTTAB | sort > $SORTEDINPUT; } || {
	echo
	echo "Error! Sorting of unique transcripts failed. Aborting. Check if"
	echo "$INPUTFILE is correct FASTA file and check if file $INPUTTAB is correct."
	echo
	exit 1
	}

# Apply the join command
join -j1 $UNIQUELIST $SORTEDINPUT > $JOINEDTS || {
	echo
	echo "Error! Joining failed. Aborting. Check files"
	echo "$UNIQUELIST and $SORTEDINPUT if they have same number of lines."
	echo
	exit 1
	}

# Convert to FASTA
echo
echo "Converting to FASTA"
sed 's/ /\t/g' $JOINEDTS > $JOINEDTABS || {
	echo
	echo "Error! Conversion of $JOINEDTS to FASTA failed."
	echo "Aborting. Check if file  $JOINEDTS is correct."
	echo
	exit 1
	}

echo
awk '{print ">"$1"\n"$2}' $JOINEDTABS > $JOINEDFA || {
	echo
	echo "Error! Conversion of $JOINEDTS to FASTA failed."
	echo "Aborting. Check if file $JOINEDTABS is correct."
	echo
	exit 1
	}
echo "Joined transcripts written in FASTA format as"
echo "  $JOINEDFA for possible later usage."
confirmgo

# Step 2: Find genome skim data (only nuclear reads), which align to the unique transcripts

echo "Step 2 of the pipeline - removal of reads of plastid origin."

# Get rid of the chloroplast and mitochondrial reads in the genome skim data

# Chloroplast reads

# Create a reference plastome index with bowtie2-build
echo
echo "Creating a reference plastome index"
echo
bowtie2-build $REFERENCECP $REFERENCECP2 || {
	echo
	echo "Error! Creating a reference plastome index with bowtie2-build failed.  Aborting."
	echo "Check if file $REFERENCECP is correct."
	echo
	exit 1
	}
echo

# Map the cpDNA reads to the reference plastome with Bowtie2
echo "Mapping cpDNA reads to the reference plastome. This may take longer time."
echo
bowtie2 -x $REFERENCECP2 -1 $INPUTFQ1 -2 $INPUTFQ2 -S $BOWTIE2CP || {
	echo
	echo "Error! Mapping cpDNA reads to the reference plastome with bowtie2 failed."
	echo "Aborting. Check if files $REFERENCECP,"
	echo "$REFERENCECP2, $INPUTFQ1"
	echo "and $INPUTFQ2 are correct."
	echo
	exit 1
	}
echo
echo "Mapping finished"

# Convert SAM to FASTQ with SAMtools
echo
echo "Converting SAM to FASTQ. This may take longer time."
samtools view -b -T $REFERENCECP $BOWTIE2CP | samtools fastq -n - -1 $FASTQNOCP.1.fq -2 $FASTQNOCP.2.fq || {
	echo
	echo "Error! Conversion of SAM to FASTQ with SAMtools failed. Aborting."
	echo "Check if files $REFERENCECP and $BOWTIE2CP are correct."
	echo
	exit 1
	}
echo

echo
echo "Removed reads saved for possible later usage as"
ls -1 $FASTQNOCP*
confirmgo

# Mitochondrial reads - optional step

if [[ -n "$REFERENCEMT" && -n "$REFERENCEMT0" ]]; then

	echo "Step 3 of the pipeline - removal of reads of mitochondrial origin (optional)."

	# Create a reference mitochondriome index with bowtie2-build
	echo
	echo "Creating a reference mitochondriome index"
	echo
	bowtie2-build $REFERENCEMT $REFERENCEMT2 || {
		echo
		echo "Error! Creating of reference mitochondriome index with bowtie2-build failed."
		echo "Aborting. Check if files $REFERENCEMT and $REFERENCEMT2 are correct."
		echo
		exit 1
		}
	echo

	# Map the mtDNA reads to the reference mitochondriome with Bowtie2
	echo "Mapping mtDNA reads to reference mitochondriome. This may take longer time."
	echo
	bowtie2 -x $REFERENCEMT2 -1 $FASTQNOCP.1.fq -2 $FASTQNOCP.2.fq -S $BOWTIE2MT || {
		echo
		echo "Error! Mapping mtDNA reads to reference mitochondriome with bowtie2 failed."
		echo "Aborting. Check if files $REFERENCEMT2,"
		echo "$FASTQNOCP.1.fq and $FASTQNOCP.2.fq are correct."
		echo
		exit 1
		}
	echo

	# Convert SAM to FASTQ with SAMtools
	echo "Converting SAM to FASTQ. This may take longer time."
	samtools view -b -T $REFERENCEMT $BOWTIE2MT | samtools fastq -n - -1 $FASTQNOMT.1.fq -2 $FASTQNOMT.2.fq || {
		echo
		echo "Error! Conversion of SAM to FASTQ with SAMtools failed. Aborting. Check if files"
		echo "$REFERENCEMT and $BOWTIE2MT are correct."
		echo
		exit 1
		}
	echo

	# Combine the paired-end reads with FLASH - with mitochondrial reads
	# Sanitize output variable for usage by FLASh (it always starts with './', thus absolute URL is not usable)
	FLASHOUT2=$(realpath --relative-to=$WORKDIR $FLASHOUT)
	echo "Step 4 of the pipeline - combination of paired-end reads."
	echo
	echo "Combining paired-end reads"
	echo
	flash -o $FLASHOUT2 -M $FLASHM $FASTQNOMT.1.fq $FASTQNOMT.2.fq || {
		echo
		echo "Error! Combining paired-end reads failed. Aborting. Check if files"
		echo "$REFERENCEMT, $FASTQNOMT.1.fq and $FASTQNOMT.2.fq are correct."
		echo
		exit 1
		}
	echo
	else
		# Combine the paired-end reads with FLASH - without mitochondrial reads
		echo "Step 4 of the pipeline - combination of paired-end reads."
		echo
		echo "Combining paired-end reads"
		echo
		flash -o $FLASHOUT2 -M $FLASHM `echo $FASTQNOCP`_1.fq `echo $FASTQNOCP`_2.fq || {
			echo
			echo "Error! Combining paired-end reads failed. Aborting. Check if files"
			echo "${REDF}$REFERENCECP${NORM}, ${REDF}`echo $FASTQNOCP`_1.fq${NORM} and ${REDF}`echo $FASTQNOCP`_2.fq${NORM} are correct."
			echo
			exit 1
			}
		fi

# Step 5 - matching of the unique transcripts and the filtered, combined genome skim reads sharing ≥85% sequence similarity

echo
echo "Step 5 of the pipeline - matching of the unique transcripts and the filtered,"
echo "  combined genome skim reads sharing ≥85% sequence similarity."

# Convert FASTQ file to FASTA
echo
echo "Converting FASTQ to FASTA. This may take longer time."
# print the first and 2nd line of every 4 lines # AWK: cat IN.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > OUT.fa
sed -n '1~4s/^@/>/p;2~4p' $FLASHOUT.extendedFrags.fastq > $FLASHOUT.extendedFrags.fa || {
	echo
	echo "Error! Conversion of FASTQ to FASTA failed. Aborting. Check if file"
	echo "$FLASHOUT.extendedFrags.fastq is correct."
	echo
	exit 1
	}
echo "Converted FASTA saved as"
echo "  $FLASHOUT.extendedFrags.fa for possible later usage."
confirmgo

# BLAT between the unique transcripts and the genome skim data
echo "BLAT between the unique transcripts and the genome skim data."
echo "This may take longer time."
blat -t=dna -q=dna -minIdentity=$BLATIDENT -out=pslx $JOINEDFA $FLASHOUT.extendedFrags.fa $BLATOUTFIN || {
	echo
	echo "Error! BLAT between the unique transcripts and the genome skim data failed."
	echo "Aborting. Check if files $JOINEDFA and $FLASHOUT.extendedFrags.fa are correct."
	echo
	exit 1
	}
echo " BLAT output saved as"
echo "  $BLATOUTFIN for possible later usage."
confirmgo

# Step 6: Assemble the obtained sequences in contigs

echo "Step 6 of the pipeline - filtering of BLAT output."

# Modification of the PSLX file is needed: remove headers (first 5 lines), select the field with the transcript (target) sequence names and the field with the query sequences, remove empty sequences
echo
echo "Modifying PSLX BLAT output for usage in Geneious"
{ sed 1,5d $BLATOUTFIN | cut -f 14,$PSLXCUT | awk '{n=split($2,a,",");for(i=1;i<=n;i++)print $1"_"NR"_"i,a[i]}' | grep "[0-9_]\+[[:blank:]][acgtuwsmkrybdhvn\-]\+" > $TAB; } || {
	echo
	echo "Error! Modifying PSLX BLAT output failed. Aborting."
	echo "Check if file $BLATOUTFIN is correct."
	echo
	exit 1
	}
echo
echo "Modified file saved as"
echo "  $TAB for possible later usage."
confirmgo

{ awk '{print $1"\t"length($2)"\t"$2}' $TAB | awk '{sum+=$2}END{print sum}'; } || {
	echo
	echo "Error! Conversion of FASTA to TAB failed. Aborting."
	echo "Check if file $TAB is correct."
	echo
	exit 1
	}

# Remove transcripts with >1000 BLAT scores (or another value selected by user; very likely that these are repetitive elements)

# Count the number of times each transcript hit a genome skim read
echo
echo "Counting number of times each transcript hit a genom skim read"
{ cut -f 1 -d "_" $TAB | sort | uniq -c | sort -n -r > $TABLIST; } || {
	echo
	echo "Error! Counting of number of times each transcript hit a genom skim read failed."
	echo "Aborting. Check if file $TAB is correct."
	echo
	exit 1
	}

# List of the transcripts with >1000 BLAT scores (or another value selected by user)
echo
echo "Listing transcripts with >$BLATSCORE BLAT scores"
{ awk '$1>'"$BLATSCORE"'' $TABLIST | awk '{print $2}' > $TABBLAT; } || {
	echo
	echo "Error! Listing of transcripts with >$BLATSCORE BLAT scores failed. Aborting."
	echo "Check if file $TABLIST is correct."
	echo
	exit 1
	}

# Make a new TAB file without these transcripts
echo
echo "Removing transcripts with >$BLATSCORE BLAT score"
grep -v -f $TABBLAT $TAB > $TABREMOVED || {
	echo
	echo "Error! Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
	echo "Check if files $TABBLAT and $TAB are correct."
	echo
	exit 1
	}
echo

{ awk '{print $1"\t"length($2)"\t"$2}' $TABREMOVED | awk '{sum+=$2}END{print sum}'; } || {
	echo
	echo "Error! Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
	echo
	exit 1
	}

# Convert this TAB file to FASTA and remove the sequences with 'n'
echo
echo "Converting TAB to FASTA and removing sequences with \"n\""
{ grep -v n $TABREMOVED | sed 's/^/>/' | sed 's/\t/\n/' > $FINALA; } || {
	echo
	echo "Error! Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
	echo "Check if file $TABREMOVED is correct."
	echo
	exit 1
	}

grep -v n $TABREMOVED | awk '{print $1"\t"length($2)}' | awk '{s+=$2;a++}END{print s}' || {
	echo
	echo "Error! Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
	echo "Check if file $TABREMOVED is correct."
	echo
	exit 1
	}

# Remove unneeded temporal files - keep only *.pslx, *.fasta and *.bam
echo
echo "Removing unneeded temporal files"
if [ -n "$REFERENCEMT0" ]; then
	rm $INPUTFILE0 $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $REFERENCECP $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $REFERENCEMT2* $REFERENCEMT $BOWTIE2M$FLASHOUT.extendedFrags.fastq $TAB $TABLIST $TABBLAT $TABREMOVED || {
		echo
		echo "Error! Removal of temporal files failed. Remove following files manually:"
		echo "  \"$INPUTFILE0\", \"$UNIQUELIST\", \"$INPUTTAB\","
		echo "  \"$SORTEDINPUT\", \"$REFERENCECP\", \"$JOINEDTS\","
		echo "  \"$JOINEDTABS\", \"$REFERENCECP2*\", \"$REFERENCEMT\","
		echo "  \"$BOWTIE2CP\", \"$REFERENCEMT2*\", \"$BOWTIE2MT\","
		echo "  \"$FLASHOUT.extendedFrags.fastq\", \"$TAB\","
		echo "  \"$TABLIST\",, \"$TABBLAT\" and \"$TABREMOVED\"."
		confirmgo
		}
	else
		rm $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $FLASHOUT.extendedFrags.fastq $TAB $TABLIST $TABBLAT $TABREMOVED || {
			echo
			echo "Error! Removal of temporal files failed. Remove following files manually:"
			echo "  \"$INPUTFILE0\", \"$UNIQUELIST\", \"$INPUTTAB\","
			echo "  \"$SORTEDINPUT\", \"$JOINEDTS\", \"$JOINEDTABS\","
			echo "  \"$REFERENCECP\", \"$REFERENCECP2*\",\"$BOWTIE2CP\","
			echo "  \"$FLASHOUT.extendedFrags.fastq\", \"$TAB\","
			echo "  \"$TABLIST\", \"$TABBLAT\" and \"$TABREMOVED\"."
			confirmgo
			}
		fi

# List kept files which user can use for another analysis
echo
echo "Following files are kept for possible later usage (see manual for details):"
echo "================================================================================"
if [ -n "$REFERENCEMT" ]; then
	echo "1) List of old names of FASTA sequences in $INPUTFILE0 and"
	echo "   renamed FASTA sequences in $INPUTFILE:"
	echo "$TRANSCRIPTOMEFASTANAMES"
	echo "2) Output of BLAT (removal of transcripts sharing ≥90% sequence similarity):"
	echo "$BLATOUT"
	echo "3) Unique transcripts in FASTA format:"
	echo "$JOINEDFA"
	echo "4) Genome skim data without cpDNA reads:"
	ls $FASTQNOCP*
	echo "5) Genome skim data without mtDNA reads:"
	ls $FASTQNOMT*
	echo "6) Combined paired-end genome skim reads:"
	echo "$FLASHOUT.extendedFrags.fa"
	echo "7) Output of BLAT (matching of the unique transcripts and the filtered,"
	echo "   combined genome skim reads sharing ≥85% sequence similarity):"
	echo "$BLATOUTFIN"
	echo "8) Matching sequences in FASTA:"
	echo "$BLATOUTFIN2"
	echo "9) Final FASTA sequences for usage in Geneious:"
	echo "$FINALA"
	else
		echo "1)  List of old names of FASTA sequences in $INPUTFILE0 and"
		echo "    renamed FASTA sequences in $INPUTFILE:"
		echo "$TRANSCRIPTOMEFASTANAMES"
		echo "2) Output of BLAT (removal of transcripts sharing ≥90% sequence similarity):"
		echo "$BLATOUT"
		echo "3) Unique transcripts in FASTA format:"
		echo "$JOINEDFA"
		echo "4) Genome skim data without cpDNA reads:"
		ls $FASTQNOCP*
		echo "5) Combined paired-end genome skim reads:"
		echo "$FLASHOUT.extendedFrags.fa"
		echo "6) Output of BLAT (matching of the unique transcripts and the filtered,"
		echo "   combined genome skim reads sharing ≥85% sequence similarity):"
		echo "$BLATOUTFIN"
		echo "7) Matching sequences in FASTA:"
		echo "$BLATOUTFIN2"
		echo "8) Final FASTA sequences for usage in Geneious:"
		echo "$FINALA"
		fi
echo "================================================================================"
confirmgo

echo "Success!"
echo
echo "================================================================================"
echo "Resulting FASTA was saved as"
echo "$FINALA"
echo "for usage in Geneious (step 7 of the pipeline)."
echo "Use this file in next step of the pipeline. See PDF manual for details."
echo "================================================================================"
confirmgo

echo "================================================================================"
echo "Run Geneious (tested with versions 6-9), see PDf manual for details"
echo
echo "Use exported files from Geneious as input for part B of the Sondovač script:"
echo "  ./sondovac_part_b.sh -h, see PDF manual for details."
echo "================================================================================"
echo
echo "Script exited successfully..."
echo

exit
