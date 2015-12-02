#!/bin/bash

# Determine script's directory
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load aliases to replace Mac OS X outdated tools by those installed by Homebrew
shopt -s expand_aliases
source $SCRIPTDIR/mac_aliases

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

echo "${REDF}This is part A of the pipeline.${NORM}"
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
OUTPUTFILENAME="output"
# Usage of transcriptome (22) or genome skim data (23, parameter "-g")
PSLXCUT=22

# Create empty variables for file names
INPUTFILE=""
INPUTFILE0=""
REFERENCECP=""
INPUTFQ1=""
INPUTFQ2=""
REFERENCEMT=""

# Parse initial arguments
while getopts "hvulrpeo:inf:c:m:t:q:a:y:s:g" START; do
  case "$START" in
    h|v)
      generaloptions
      echo
      echo -e "\tIf options ${BOLD}-f, -c, -m, -t${NORM} and/or ${BOLD}-q${NORM} are used and the script is running"
      echo -e "\t  in interactive mode, those values will be used as defaults, but may"
      echo -e "\t  later be overwritten."
      echo
      echo -e "\tOptions required for running in non-interactive mode:"
      echo -e "\t${REDF}-f${NORM}\t${CYAF}Transcriptome input file${NORM} in FASTA format."
      echo -e "\t${REDF}-c${NORM}\t${CYAF}Plastome reference sequence${NORM} input file in FASTA format."
      echo -e "\t${REDF}-m${NORM}\t${CYAF}Mitochondriome reference sequence${NORM} input file in FASTA format."
      echo -e "\t\t  This file is optional. In interactive mode you will be every"
      echo -e "\t\t  time asked if you wish to use it."
      echo -e "\t${REDF}-t${NORM}\t${CYAF}Paired-end genome skim input file${NORM} in FASTQ format (first file)."
      echo -e "\t${REDF}-q${NORM}\t${CYAF}Paired-end genome skim input file${NORM} in FASTQ format (second file)."
      echo
      echo -e "\tOther optional arguments (if not provided, default values are used):"
      echo -e "\t${REDF}-a${NORM}\t${CYAF}Maximum overlap length expected in approximately 90% of read pairs${NORM} (parameter \"-M\" of"
      echo -e "\t\t  FLASH, see its manual for details)."
      echo -e "\t\tDefault value: 65 (integer ranging from 10 to 300)"
      echo -e "\t${REDF}-y${NORM}\t${CYAF}Sequence similarity between unique transcripts and the filtered,"
      echo -e "\t\t  combined genome skim reads${NORM} (parameter \"-minIdentity\" of BLAT,"
      echo -e "\t\t  see its manual for details)."
      echo -e "\t\tDefault value: 85 (integer ranging from 70 to 100; the default value of 85% minimum sequence similarity suggests gene orthology)"
      echo -e "\t${REDF}-s${NORM}\t${CYAF}Number of BLAT hits per transcript${NORM} when matching unique"
      echo -e "\t\t  transcripts and the filtered, combined genome skim reads."
      echo -e "\t\tDefault value: 1000 (integer ranging from 100 to 10000)"
      echo -e "\t${REDF}-g${NORM}\t${CYAF}Use genome skim sequences instead of transcripts${NORM} for making the"
      echo -e "\t\t  probes. Default is usage of transcripts (no parameter)."
      echo -e "\t${BOLD}WARNING!${NORM} If parameters ${BOLD}-a, -y${NORM}, ${BOLD}-s${NORM} or ${BOLD}-g${NORM} are not provided, default values"
      echo -e "\t\t are taken and it is not possible to change them later (not even"
      echo -e "\t\t in interactive mode)."
      echo
      echo "Examples:"
      echo "Basic and the most simple usage:"
      echo "${REDF}$0${NORM} ${CYAF}-i${NORM}"
      echo "Specify some of required input files, otherwise run interactively:"
      echo "${REDF}$0${NORM} ${CYAF}-i -f${NORM} input.fa ${CYAF}-t${NORM} reads1.fastq ${CYAF}-q${NORM} reads2.fastq"
      echo "Running in non-interactive automated way:"
      echo "${REDF}$0${NORM} ${CYAF}-n -f${NORM} input.fa ${CYAF}-c${NORM} referencecp.fa ${CYAF}-m${NORM} referencemt.fa ${CYAF}-t${NORM} reads1.fastq ${CYAF}-q${NORM} reads2.fastq"
      echo "Modify parameter ${BOLD}-a${NORM}, otherwise run interactively:"
      echo "${REDF}$0${NORM} ${CYAF}-i -a${NORM} 300"
      echo "Run in non-interactive mode (parameter ${BOLD}-n${NORM}) - in such case user must specify all"
      echo "  required input files (parameters ${BOLD}-f${NORM}, ${BOLD}-c${NORM}, ${BOLD}-m${NORM}, ${BOLD}-t${NORM} and ${BOLD}-q${NORM}). Moreover, parameter"
      echo "  ${BOLD}-y${NORM} is modified:"
      echo "${REDF}$0${NORM} ${CYAF}-n -f${NORM} input.fa ${CYAF}-c${NORM} referencecp.fa ${CYAF}-m${NORM} referencemt.fa ${CYAF}-t${NORM} reads1.fastq ${CYAF}-q${NORM} reads2.fastq ${CYAF}-y${NORM} 90"
      echo "Modifying parameter ${BOLD}-s${NORM}. Note interactive mode ${BOLD}-i${NORM} is implicit and does not need"
      echo "  to be specified explicitly:"
      echo "${REDF}$0${NORM} ${CYAF}-s${NORM} 950"
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
      INPUTFILE0=$OPTARG
      echo "Input file: ${REDF}$INPUTFILE0${NORM}"
      ;;
    c)
      REFERENCECP=$OPTARG
      echo "Plastom reference: ${REDF}$REFERENCECP${NORM}"
      ;;
    m)
      REFERENCEMT=$OPTARG
      echo "Chondrion reference: ${REDF}$REFERENCEMT${NORM}"
      ;;
    t)
      INPUTFQ1=$OPTARG
      echo "FASTQ reads 1: ${REDF}$INPUTFQ1${NORM}"
      ;;
    q)
      INPUTFQ2=$OPTARG
      echo "FASTQ reads 2: ${REDF}$INPUTFQ2${NORM}"
      ;;
    a)
      FLASHM=$OPTARG
      # Check if provided value makes sense
      case "$FLASHM" in
	125) FLASHM=125;;
	150) FLASHM=150;;
	250) FLASHM=250;;
	300) FLASHM=300;;
	*) echo
	  echo "${REDF}${BOLD}Error!${NORM} For parameter \"-a\" you did not provide any of the allowed values 125, 150,"
	  echo "   250, or 300!"
	  echo
	  exit 1
	  ;;
	esac
      echo "Read length of paired-end reads: ${REDF}$FLASHM${NORM}"
      ;;
    y)
      BLATIDENT=$OPTARG
      # Check if provided value makes sense
      if [[ "$BLATIDENT" =~ ^[0-9]+$ ]] && [ "$BLATIDENT" -ge 70 -a "$BLATIDENT" -le 100 ]; then
	echo "BLAT score for identity between unique transcripts and genome skim data: ${REDF}$BLATIDENT${NORM}"
	else
	  echo
	  echo "${REDF}${BOLD}Error!${NORM} For parameter \"-y\" you did not provide an integer of range from 70 to 100!"
	  echo
	  exit 1
	fi
      ;;
    s)
      BLATSCORE=$OPTARG
      # Check if provided value makes sense
      if [[ "$BLATSCORE" =~ ^[0-9]+$ ]] && [ "$BLATSCORE" -ge 100 -a "$BLATSCORE" -le 10000 ]; then
	echo "BLAT score: ${REDF}$BLATSCORE${NORM}"
	else
	  echo
	  echo "${REDF}${BOLD}Error!${NORM} For parameter \"-s\" you did not provide an integer ranging from 100 to 10000!"
	  echo
	  exit 1
	fi
      ;;
    g)
      PSLXCUT=23
      echo "${CYAF}The script will use genome skim data instead of transcriptome.${NORM}"
      ;;
    ?)
      echo
      echo "Invalid option(s)!"
      echo "See \"${REDF}$0 -h${NORM}\" for usage options."
      echo
      exit 1
      ;;
    esac
  done

# Check if user didn't use together -n and -i
checkmodef

# Ensure user reads introductory information
confirmgo

# NOTE: warn user this is not yet for daily usage
devrelease

# Check operating system
oscheck

# Set variables for working directory and PATH
workdirpath

# Check availability of all needed binaries

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
BOWTIE2V="2.2.6"
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
  echo "\"${REDF}bowtie2${NORM}\", \"${REDF}bowtie-align-l${NORM}\", \"${REDF}bowtie-align-s${NORM}\", \"${REDF}bowtie2-build${NORM}\", \"${REDF}bowtie2-build-l${NORM}\"" &&
  echo "   and \"${REDF}bowtie2-build-s${NORM}\" are available. ${GREF}OK.${NORM}"
  } || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
    echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
    echo "  download bowtie2-$BOWTIE2V-source.zip, compile it and ensure it is in PATH."
    echo
    exit 1
    }
  }

# Check if bowtie2 is available
{ command -v bowtie2 >/dev/null 2>&1 &&
  command -v bowtie2-align-l >/dev/null 2>&1 &&
  command -v bowtie2-align-s >/dev/null 2>&1 &&
  command -v bowtie2-build >/dev/null 2>&1 &&
  command -v bowtie2-build-l >/dev/null 2>&1 &&
  command -v bowtie2-build-s >/dev/null 2>&1 &&
  echo "\"${REDF}bowtie2${NORM}\", \"${REDF}bowtie-align-l${NORM}\", \"${REDF}bowtie-align-s${NORM}\", \"${REDF}bowtie2-build${NORM}\", \"${REDF}bowtie2-build-l${NORM}\"" &&
  echo "  and \"${REDF}bowtie2-build-s${NORM}\" are available. ${GREF}OK.${NORM}"
  } || {
  echo
  echo "\"${CYAF}bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\""
  echo "  and \"bowtie2-build-s\" are required but not installed or available in PATH.${NORM}"
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile Bowtie2-$BOWTIE2V from source available together with this script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to compile Bowtie2-$BOWTIE2V from source code${NORM} downloaded from"
    echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
    echo "Type \"${REDF}D${NORM}\" ${CYAF}to download Bowtie2-$BOWTIE2V binary${NORM} from"
    echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM} for your OS."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy Bowtie2-2.$BOWTIE2V binary available together with the script${NORM}"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
    echo "  See \"${REDF}brew info homebrew/science/bowtie2${NORM}\" for more details."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to install"
    echo "  Bowtie2 yourselves."
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
	    echo "Downloading Bowtie2 binaries for $OS $OSB" &&
	    $DOWNLOADER bowtie2-$BOWTIE2V-linux-x86_64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-linux-x86_64.zip &&
	    unzip -nq bowtie2-$BOWTIE2V-linux-x86_64.zip &&
	    cp bowtie2-$BOWTIE2V/bowtie2* $BIN/ &&
	    echo "\"${REDF}bowtie2${NORM}\", \"${REDF}bowtie-align-l${NORM}\", \"${REDF}bowtie-align-s${NORM}\", \"${REDF}bowtie2-build${NORM}\", \"${REDF}bowtie2-build-l${NORM}\"" &&
	    echo "  and \"${REDF}bowtie2-build-s${NORM}\" are available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to"
	      echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
	      echo "  and download and unpack bowtie2-$BOWTIE2V-linux-x86_64.zip manually."
	      echo
	      exit 1
	      }
	  elif [ "$OS" == "Mac" ]; then
	    {
	    echo "Downloading Bowtie2 binaries for $OS $OSB" &&
	    $DOWNLOADER bowtie2-$BOWTIE2V-macos-x86_64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-macos-x86_64.zip &&
	    unzip -nq bowtie2-$BOWTIE2V-macos-x86_64.zip &&
	    cp bowtie2-$BOWTIE2V/bowtie2* $BIN/ &&
	    echo "\"${REDF}bowtie2${NORM}\", \"${REDF}bowtie-align-l${NORM}\", \"${REDF}bowtie-align-s${NORM}\", \"${REDF}bowtie2-build${NORM}\", \"${REDF}bowtie2-build-l${NORM}\"" &&
	    echo "  and \"${REDF}bowtie2-build-s${NORM}\" are available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to"
	      echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
	      echo "  and download and unpack bowtie2-$BOWTIE2V-macos-x86_64.zip manually."
	      echo
	      exit 1
	      }
	  elif [ "$OS" == "Windows" ]; then
	    {
	    echo "Downloading Bowtie2 binaries for $OS $OSB" &&
	    $DOWNLOADER bowtie2-$BOWTIE2V-mingw-win64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-mingw-win64.zip &&
	    unzip -nq bowtie2-$BOWTIE2V-mingw-win64.zip &&
	    cp bowtie2-$BOWTIE2V/bowtie2* $BIN/ &&
	    echo "\"${REDF}bowtie2${NORM}\", \"${REDF}bowtie-align-l${NORM}\", \"${REDF}bowtie-align-s${NORM}\", \"${REDF}bowtie2-build${NORM}\", \"${REDF}bowtie2-build-l${NORM}\""
	    echo "  and \"${REDF}bowtie2-build-s${NORM}\" are available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to"
	      echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
	      echo "  and download and unpack bowtie2-$BOWTIE2V-mingw-win64.zip manually."
	      echo
	      exit 1
	      }
	  else
	    echo "Unknown OS or OS without Bowtie2 binary available."
	    compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
	  fi
	  break
	  ;;
	S|s)
	  echo
	  downloaderselector
	  checktools unzip
	  echo
	  echo "Downloading Bowtie2 source code"
	  $DOWNLOADER bowtie2-$BOWTIE2V-source.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$BOWTIE2V/bowtie2-$BOWTIE2V-source.zip || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to"
	    echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
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
	      echo "Copying Bowtie binaries"
	      cp -pr $SCRIPTDIR/pkgs/macosx/bin/bowtie2* $BIN/
	      mkdir -p $WORKDIR/bin/share
	      cp -pr $SCRIPTDIR/pkgs/macosx/share/man $WORKDIR/bin/share/
	      break
	      ;;
	    Linux)
	      echo "Copying Bowtie binaries"
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
	    { echo "Installing ${REDF}Bowtie2${NORM} using Homebrew" &&
	    brew install homebrew/science/bowtie2 &&
	    echo "\"${REDF}Bowtie2${NORM}\" is available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of${NORM} \"${REDF}Bowtie2${NORM}\" ${CYAF}failed.${NORM} Please, do it manually. For details see"
	      echo "\"${REDF}brew info homebrew/science/bowtie2${NORM}\" and \"${REDF}brew help${NORM}\"."
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
	  echo "Please, go to ${REDF}http://bowtie-bio.sourceforge.net/bowtie2/index.shtml${NORM} and install"
	  echo " latest Bowtie2 and ensure it is in PATH."
	  exit 2
	  ;;
	*) echo "${CYAF}Wrong option.${NORM} Use ${REDF}C${NORM}, ${REDF}D${NORM}, ${REDF}S${NORM}, ${REDF}B${NORM}, ${REDF}H${NORM} or ${REDF}M${NORM}." && read BOWTIE;;
      esac
    done
  else
    exit 1
  fi
  }

function compilesamtools {
  cd $SCRIPTDIR/src/samtools-1.2 &&
  echo &&
  checktools make &&
  echo &&
  echo "Compiling samtools..." &&
  make -s &&
  make -s prefix=$WORKDIR/bin install &&
  cd $WORKDIR &&
  echo "\"${REDF}samtools${NORM}\" is available. ${GREF}OK.${NORM}"
  }

# Check if samtools is available
{ command -v samtools >/dev/null 2>&1 && echo "\"${REDF}samtools${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "\"${REDF}samtools${NORM}\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile SAMtools-1.2 from source available together with this script.${NORM}"
    echo "  Makefile was modified not to require GNU ncurses library."
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download latest developmental SAMtools source${NORM} from"
    echo "  ${REDF}https://github.com/samtools/samtools/${NORM} and compile it. Compilation requires GNU"
    echo "  ncurses library and is recommended only for advanced users. If compilation"
    echo "  fails, check SAMtools' ${REDF}INSTALL${NORM} file for details and adjust its Makefile."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy SAMtools-1.2 binary${NORM} available together with the script"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
    echo "  See \"${REDF}brew info homebrew/science/samtools${NORM}\" for more details."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install SAMtools yourselves."
    read SAMTOOLS
    while :
    do
      case "$SAMTOOLS" in
	C|c)
	  compilesamtools || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to ${REDF}http://www.htslib.org/download/${NORM}"
	    echo "  download samtools-1.2, compile it and ensure it is in PATH."
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
	  echo "Downloading SAMtools sources..." &&
	  $DOWNLOADER develop.zip https://github.com/samtools/samtools/archive/develop.zip &&
	  unzip -nq develop.zip &&
	  cd samtools-develop &&
	  $DOWNLOADER develop.zip https://github.com/samtools/htslib/archive/develop.zip &&
	  unzip -nq develop.zip &&
	  echo "Compiling SAMtools" &&
	  make -s HTSDIR=`pwd`/htslib-develop &&
	  make -s prefix=$WORKDIR/bin install &&
	  cd $WORKDIR &&
	  echo "\"${REDF}samtools${NORM}\" is available. ${GREF}OK.${NORM}"
	  } || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to ${REDF}http://www.htslib.org/download/${NORM} for"
	    echo "  latest stable version of SAMtools or ${REDF}https://github.com/samtools/${NORM} for latest"
	    echo "  developmental version, download samtools, compile it and ensure it is in PATH."
	    echo
	    exit 1
	    }
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying SAMtools binaries"
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
	      echo "Copying SAMtools binaries"
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
	  { echo "Installing ${REDF}SAMtools using Homebrew" &&
	  brew install homebrew/science/samtools &&
	  echo "\"${REDF}SAMtools${NORM}\" is available. ${GREF}OK.${NORM}"
	  } || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of${NORM} \"${REDF}SAMtools${NORM}\" ${CYAF}failed.${NORM} Please, do it manually. For details see"
	    echo "\"${REDF}brew info homebrew/science/samtools${NORM}\" and \"${REDF}brew help${NORM}\"."
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
	  echo "Please, go to ${REDF}http://www.htslib.org/${NORM} and install SAMtools 1.2 and ensure it is"
	  echo "  in PATH."
	  exit 2
	  ;;
	*) echo "${CYAF}Wrong option.${NORM} Use ${REDF}C${NORM}, ${REDF}S${NORM}, ${REDF}B${NORM}, ${REDF}H${NORM} or ${REDF}M${NORM}." && read SAMTOOLS;;
      esac
    done
  else
    exit 1
  fi
  }

function compilebam2fastq {
  echo &&
  checktools make &&
  checktools g++ &&
  echo &&
  echo "Compiling bam2fastq..." &&
  cd $1 &&
  make -s &&
  cp bam2fastq $BIN/ &&
  cd $WORKDIR &&
  echo "\"${REDF}bam2fastq${NORM}\" is available. ${GREF}OK.${NORM}"
  }

# Check if bam2fastq is available
{ command -v bam2fastq >/dev/null 2>&1 && echo "\"${REDF}bam2fastq${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "\"${REDF}bam2fastq${NORM}\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile bam2fastq from source available together with this script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download bam2fastq source${NORM} from"
    echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM} and compile it."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy bam2fastq binary available together with the script${NORM}"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install bam2fastq yourselves."
    read bam2fastq
    while :
    do
      case "$bam2fastq" in
	C|c)
	  compilebam2fastq $SCRIPTDIR/src/bam2fastq-1.1.0 || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
	    echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}, download bam2fastq,"
	    echo "  compile it manually and ensure it is in PATH."
	    echo
	    exit 1
	    }
	  break
	  ;;
	S|s)
	  {
	  echo &&
	  downloaderselector &&
	  checktools tar &&
	  checktools gunzip &&
	  checktools make &&
	  echo &&
	  echo "Downloading bam2fastq source code..." &&
	  $DOWNLOADER bam2fastq-1.1.0.tgz http://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz &&
	  tar xzvf bam2fastq-1.1.0.tgz &&
	  compilebam2fastq bam2fastq-1.1.0/
	  } || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
	    echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}, download bam2fastq,"
	    echo "  compile it manually and ensure it is in PATH."
	    echo
	    exit 1
	    }
		break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying bam2fastq binary"
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/bam2fastq $BIN/
	      break
	      ;;
	    Linux)
	      echo "Copying bam2fastq binary"
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/bam2fastq $BIN/
	      break
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB."
	      echo
	      compilebam2fastq $SCRIPTDIR/src/bam2fastq-1.1.0 || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
		echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}, download bam2fastq,"
		echo "  compile it manually and ensure it is in PATH."
		echo
		exit 1
		}
	  break
	      ;;
	  esac
	  break
	  ;;
	M|m)
	  echo "Please, go to ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}"
	  echo "  download bam2fastq, compile it manually and ensure it is in PATH."
	  exit 2
	  ;;
	*) echo "${CYAF}Wrong option.${NORM} Use ${REDF}C${NORM}, ${REDF}S${NORM}, ${REDF}B${NORM} or ${REDF}M${NORM}." && read bam2fastq;;
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
  echo "Compiling FLASH from source code..." &&
  echo &&
  cd $1 &&
  make -s &&
  cp flash $BIN/ &&
  cd $WORKDIR &&
  echo &&
  echo "\"${REDF}flash${NORM}\" is available. ${GREF}OK.${NORM}"
  } || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
    echo "  ${REDF}http://sourceforge.net/projects/flashpage/files/${NORM} download latest"
    echo "  FLASH-*.tar.gz, compile it manually and ensure it is in PATH."
    echo
    exit 1
    }
  }

# Check if FLASH is available
{ command -v flash >/dev/null 2>&1 && echo "\"${REDF}flash${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "${CYAF}FLASH is required but not installed or available in PATH.${NORM}"
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile FLASH from source available together with this script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download FLASH source${NORM} from"
    echo "  ${REDF}http://sourceforge.net/projects/flashpage/${NORM} and compile it."
    echo "Type \"${REDF}D${NORM}\" ${CYAF}to download FLASH binary${NORM} from"
    echo "  ${REDF}http://sourceforge.net/projects/flashpage/${NORM} (available only for Windows)."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy FLASH 1.2.11 binary available together with the script${NORM}"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
    echo "  See \"${REDF}brew info homebrew/science/flash${NORM}\" for more details."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install FLASH yourselves."
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
	  echo "Downloading FLASH source code"
	  $DOWNLOADER FLASH-1.2.11.tar.gz http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to ${REDF}http://sourceforge.net/projects/flashpage/files/${NORM}"
	    echo "  download FLASH-*.tar.gz, compile it manually and ensure it is in PATH."
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
	      echo "Downloading FLASH for $OS"
	      $DOWNLOADER FLASH-1.2.11-windows-bin.zip http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11-windows-bin.zip || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to ${REDF}http://sourceforge.net/projects/flashpage/files/${NORM}"
		echo "  download FLASH-*windows-bin.zip, unpack it and ensure it is in PATH."
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
	      echo "Copying FLASH binary"
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/flash $BIN/
	      break
	      ;;
	    Linux)
	      echo "Copying FLASH binary"
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
	    { echo "Installing FLASH using Homebrew" &&
	    brew install homebrew/science/flash &&
	    echo "\"${REDF}FLASH${NORM}\" is available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of${NORM} \"${REDF}FLASH${NORM}\" ${CYAF}failed.${NORM} Please, do it manually. For details see"
	      echo "\"${REDF}brew info homebrew/science/flash${NORM}\" and \"${REDF}brew help${NORM}\"."
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
	  echo "Please, go to ${REDF}http://ccb.jhu.edu/software/FLASH/${NORM} and install FLASH and ensure it"
	  echo "  is in PATH"
	  exit 2
	  ;;
	*) echo "${CYAF}Wrong option.${NORM} Use ${REDF}C${NORM}, ${REDF}S${NORM}, ${REDF}D${NORM}, ${REDF}B${NORM}, ${REDF}H${NORM} or ${REDF}M${NORM}." && read FLASH;;
      esac
    done
  else
    exit 1
  fi
  }

# Function to compile FASTX-Toolkit
function compilefastx {
  {
  echo &&
  checktools make &&
  checktools g++ &&
  checktools realpath &&
  checktools pkg-config &&
  echo &&
  echo "Compiling FASTX-Toolkit from source code..." &&
  echo &&
  cd $1 &&
  ./configure --prefix=$WORKDIR/bin &&
  make -s &&
  make -s check &&
  make -s install &&
  make -s installcheck &&
  export PKG_CONFIG_PATH=`realpath $WORKDIR/bin/*/pkgconfig` &&
  cd $2 &&
  ./configure --prefix=$WORKDIR/bin &&
  make -s &&
  make -s check &&
  make -s install &&
  make -s installcheck &&
  cd $WORKDIR &&
  echo &&
  echo "\"${REDF}fastq_to_fasta${NORM}\" is available. ${GREF}OK.${NORM}"
  } || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
    echo "  ${REDF}http://hannonlab.cshl.edu/fastx_toolkit/download.html${NORM}, download latest"
    echo "  libgtextutils-*.tar.gz and fastx_toolkit-*.tar.bz2, compile them and ensure"
    echo "  they are in PATH."
    echo
    exit 1
    }
  }

# Check if fastq_to_fasta is available
{ command -v fastq_to_fasta >/dev/null 2>&1 && echo "\"${REDF}fastq_to_fasta${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "\"${REDF}fastq_to_fasta${NORM}\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"C\" to compile FASTX-Toolkit from source available together with this"
    echo "  script."
    echo "Type \"S\" to download FASTX-Toolkit source from"
    echo "  http://hannonlab.cshl.edu/fastx_toolkit/ and compile it."
    echo "Type \"B\" to copy FASTX-Toolkit 0.0.14 binary available together with the"
    echo "  script (recommended, available for Linux and Mac OS X)."
    echo "Type \"H\" for installation using Homebrew (only for Mac OS X, recommended)."
    echo "  See \"brew info homebrew/science/fastx_toolkit\" for more details."
    echo "Type \"M\" for manual installation - script will exit and you will have to"
    echo "  install FASTX-Toolkit yourselves."
    read FASTX
    while :
    do
      case "$FASTX" in
	C|c)
	  compilefastx $SCRIPTDIR/src/libgtextutils-0.7/ $SCRIPTDIR/src/fastx_toolkit-0.0.14/
	  break
	  ;;
	S|s)
	  echo
	  downloaderselector
	  checktools tar
	  checktools gunzip
	  checktools bunzip2 || checktools bunzip
	  echo
	  echo "Downloading FASTX source code"
	  echo
	  { $DOWNLOADER libgtextutils-0.7.tar.gz https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz &&
	  $DOWNLOADER fastx_toolkit-0.0.14.tar.bz2 https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2; } || {
	    echo
	    echo "Download failed. Please, go to"
	    echo "  https://github.com/agordon/fastx_toolkit/releases/download/ and download"
	    echo "  libgtextutils and fastx_toolkit, install them and ensure they are in PATH."
	    echo
	    exit 1
	    }
	  tar xzf libgtextutils-0.7.tar.gz
	  tar xjvf fastx_toolkit-0.0.14.tar.bz2
	  compilefastx $WORKDIR/libgtextutils-0.7/ $WORKDIR/fastx_toolkit-0.0.14/
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying FASTX binaries"
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/fasta_* $BIN/
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/fastq_* $BIN/
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/fastx_* $BIN/
	      mkdir -p $WORKDIR/bin/include $WORKDIR/bin/lib $WORKDIR/bin/share
	      cp -pr $SCRIPTDIR/pkgs/macosx/include/gtextutils $WORKDIR/bin/include/
	      cp -pr $SCRIPTDIR/pkgs/macosx/lib/pkgconfig $WORKDIR/bin/lib/
	      cp -p $SCRIPTDIR/pkgs/macosx/lib/libgtextutils* $WORKDIR/bin/lib/
	      cp -pr $SCRIPTDIR/pkgs/macosx/share/aclocal $WORKDIR/bin/share/
	      ;;
	    Linux)
	      echo "Copying FASTX binaries"
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/fasta_* $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/fastq_* $BIN/
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/fastx_* $BIN/
	      mkdir -p $WORKDIR/bin/include $WORKDIR/bin/lib64 $WORKDIR/bin/share
	      cp -pr $SCRIPTDIR/pkgs/linux64b/include/gtextutils $WORKDIR/bin/include/
	      cp -pr $SCRIPTDIR/pkgs/linux64b/lib64/pkgconfig $WORKDIR/bin/lib64/
	      cp -p $SCRIPTDIR/pkgs/linux64b/lib64/libgtextutils* $WORKDIR/bin/lib64/
	      cp -pr $SCRIPTDIR/pkgs/linux64b/share/aclocal $WORKDIR/bin/share/
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB."
	      echo
	      compilefastx $SCRIPTDIR/src/libgtextutils-0.7/ $SCRIPTDIR/src/fastx_toolkit-0.0.14/
	      ;;
	  esac
	  break
	  ;;
	H|h)
	  if [ "$OS" == "Mac" ]; then			
	    { echo "Installing FASTX Toolkit using Homebrew" &&
	    brew install homebrew/science/fastx_toolkit &&
	    echo "\"FASTX Toolkit\" is available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "Installation of \"FASTX Toolkit\" failed. Please, do it manually. For details see"
	      echo "\"brew info homebrew/science/fastx_toolkit\" and \"brew help\"."
	      echo
	      exit 1
	      }
	    else
	      echo "This is not Mac OS X. Going to compile..."
	      compilefastx $SCRIPTDIR/src/libgtextutils-0.7/ $SCRIPTDIR/src/fastx_toolkit-0.0.14/
	    fi
	  break
	  ;;
	M|m)
	  echo "Please, go to http://hannonlab.cshl.edu/fastx_toolkit/download.html download"
	  echo "  libgtextutils-*.tar.gz and fastx_toolkit-*.tar.bz2, compile them and ensure"
	  echo "  they are in PATH"
	  exit 2
	  ;;
	*) echo "Wrong option. Use C, S, B, H or M." && read FASTX;;
      esac
    done
  else
    exit 1
  fi
  }

# Input files
CHECKFILEREADOUT=""

# Input data, transcriptomic data in FASTA format
readinputfile -f "transcriptome input file in FASTA format" $INPUTFILE0
INPUTFILE0=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, plastom reference in FASTA format
readinputfile -c "plastome reference sequence input file in FASTA format" $REFERENCECP
REFERENCECP=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, HybSeq reads in FASTQ, file 1
readinputfile -m "paired-end genome skim input file in FASTQ format (first file)" $INPUTFQ1
INPUTFQ1=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, HybSeq reads in FASTQ, file 2
readinputfile -t "paired-end genome skim input file in FASTQ format (second file)" $INPUTFQ2
INPUTFQ2=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

# Input data, mitochondrion reference in FASTA
if [ "$STARTINI" == "I" ]; then
  echo
  echo "Would you like to use mitochondriome reference sequence input file in FASTA"
  echo "  format? (Yes/No)"
  read MTINPUTQ
  while :
  do
    case $MTINPUTQ in
      Y|y|Yes|yes|YES)
	readinputfile -q "mitochondriome reference sequence input file in FASTA format" $REFERENCEMT
	REFERENCEMT=$CHECKFILEREADOUT
	CHECKFILEREADOUT=""
	break
	;;
      N|n|No|no|NO)
	echo
	echo "OK, we will not use mitochondriome reference sequence. Continuing."
	REFERENCEMT=""
	break
	;;
      *) echo "Wrong option. Use Y or N." && read $MTINPUTQ;;
    esac
  done
  echo
  fi

if [ -z "$REFERENCEMT" ]; then
  echo
  echo "${BOLD}Warning!${NORM} There is no mitochondriome reference sequence."
  echo
  fi

# Input file in FASTA format
echo "Input file: $INPUTFILE0"
INPUTFILE="${INPUTFILE0%.*}_renamed.fasta"
# List of old and new names of the transcriptome FASTA sequences
TRANSCRIPTOMEFASTANAMES="${INPUTFILE0%.*}_old_and_new_names.tsv"
# Output of BLAT (removal of transcripts sharing ≥90% sequence similarity)
BLATOUT="${OUTPUTFILENAME%.*}_blat_unique_transcripts.psl"
# List of unique transcripts - temporary file - will be deleted
UNIQUELIST="${OUTPUTFILENAME%.*}_trans-trans_unique_transcripts_sorted.txt"
# Input file converted into TXT - temporary file - will be deleted
INPUTTAB="${OUTPUTFILENAME%.*}.txt"
# Sorted input file in TXT - temporary file - will be deleted
SORTEDINPUT="${OUTPUTFILENAME%.*}_sorted.txt"
# Joined unique transcriptomes - temporary file - will be deleted
JOINEDTS="${OUTPUTFILENAME%.*}_unique_transcriptomes_trans-trans_plus_sequence.txt"
# Joined unique transcriptomes in tabular format - temporary file - will be deleted
JOINEDTABS="${OUTPUTFILENAME%.*}_tabs.txt"
# Joined unique transcriptomes in FASTA format
JOINEDFA="${OUTPUTFILENAME%.*}_unique_transcripts.fasta"
# Input - reference genome - cpDNA
echo "Input file: $REFERENCECP"
# Reference genome - plastome index - temporary file - will be deleted
REFERENCECP2="${OUTPUTFILENAME%.*}.cp"
# Input cpDNA reads in FASTQ
echo "Input file: $INPUTFQ1"
echo "Input file: $INPUTFQ2"
# cpDNA reads mapped to reference - temporary file - will be deleted
BOWTIE2CP="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads.sam"
# SAM converted into BAM (removal of reads of plastid origin)
CPBAM="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads.bam"
# Genome skim data without cpDNA reads
FASTQNOCP="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads"
# Input - reference genome - mtDNA
if [ -n "$REFERENCEMT" ]; then
  echo "Input file: $REFERENCEMT"
  fi
# Reference genome - mitochondrion index - temporary file - will be deleted
REFERENCEMT2="${OUTPUTFILENAME%.*}.mt"
# mtDNA reads mapped to reference - temporary file - will be deleted
BOWTIE2MT="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_no_mt_reads.sam"
# SAM converted into BAM (removal of reads of mitochondrial origin)
MTBAM="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_no_mt_reads.bam"
# Genome skim data without mtDNA reads
FASTQNOMT="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_no_mt_reads"
# Combined paired-end genome skim reads
FLASHOUT="${OUTPUTFILENAME%.*}_combined_reads_no_cp_no_mt_reads"
# Output of BLAT (matching of the unique transcripts and the filtered, combined genome skim reads sharing ≥85% sequence similarity)
BLATOUTFIN="${OUTPUTFILENAME%.*}_blat_unique_transcripts_versus_genome_skim_data.pslx"
# Matching sequences in FASTA
BLATOUTFIN2="${OUTPUTFILENAME%.*}_blat_unique_transcripts_versus_genome_skim_data.fasta"
# FASTA converted into TABe - temporary file - will be deleted
TAB="${OUTPUTFILENAME%.*}_final.tab"
# Number of times each transcript hit a genome skim read - will be deleted
TABLIST="${OUTPUTFILENAME%.*}_transcript_hits.txt"
# Listed transcripts with >1000 BLAT score - will be deleted
TABBLAT="${OUTPUTFILENAME%.*}_1k_transcripts"
# Transcribds without >1000 BLAT score - will be deleted
TABREMOVED="${OUTPUTFILENAME%.*}_1k_transcripts-removed.tab"
# Final FASTA sequences for usage in Geneious
FINALA="${OUTPUTFILENAME%.*}_blat_unique_transcripts_versus_genome_skim_data-no_missing_fin.fsa"

# Checking if transcriptome input file has required naming scheme - only unique numbers

function renamefastanames {
  checktools paste
  echo "Input file $1"
  echo "  contains some non-allowed names. New file with correct names"
  echo "  (only increasing unique numbers) will be created."
  echo "Depending on size of your transcriptome file it can take longer time."
  while read LINE; do
    N=$((++N)) &&
    echo $LINE | sed -e 's/^>.*$/>'$N'/'
    done < $1 > $2 || {
      echo "${REDF}${BOLD}Error!${NORM} Renaming failed. Aborting. Ensure $1"
      echo "  has as names of sequences only unique numbers (nothing else)."
      echo
      exit 1
      }
  # Giving to user list of old and new names
  grep '^>' $1 > transcript_fasta_names_old
  grep '^>' $2 > transcript_fasta_names_new
  echo -e "Old FASTA names\tNew FASTA names" > $TRANSCRIPTOMEFASTANAMES
  paste transcript_fasta_names_old transcript_fasta_names_new >> $TRANSCRIPTOMEFASTANAMES
  sed -i 's/>//g' $TRANSCRIPTOMEFASTANAMES
  rm transcript_fasta_names_old transcript_fasta_names_new
  echo
  echo "File $TRANSCRIPTOMEFASTANAMES contains two columns:"
  echo "  1) Old names in original $1 and"
  echo "  2) New names in required format as in $2"
  echo "  The sequences remain intact. This might be needed to trace back some sequences."
  confirmgo
  RENAMEDFASTASEQENCES="YES"
  }

echo
echo "Checking if transcriptome input file $INPUTFILE0"
echo "  has required naming scheme - only unique numbers."
if grep '>' $INPUTFILE0 | grep -q -v '^[> ]\+[[:digit:]]\+$'; then
  echo
  renamefastanames $INPUTFILE0 $INPUTFILE
  elif  grep '>' $INPUTFILE0 | grep '^[> ]\+[[:digit:]]\+$' | uniq -c | grep -q -v '^[ ]*1'; then
    echo
    renamefastanames $INPUTFILE0 $INPUTFILE
    else
      sed -i 's/ *//g' $INPUTFILE0
      echo
      echo "Input file $INPUTFILE0"
      echo " has correct names. Continuing."
      INPUTFILE=$INPUTFILE0
      RENAMEDFASTASEQENCES="NO"
    fi

# Step 1: Obtain unique transcripts.

echo
echo "Step 1 of the pipeline - removal of transcripts sharing ≥90% sequence similarity."

# BLAT between the transcriptome itself
echo
echo "Launching BLAT between the transcriptome itself"
blat -t=dna -q=dna -noHead -out=psl $INPUTFILE $INPUTFILE $BLATOUT || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} BLAT failed. Aborting. Check if $INPUTFILE is correct."
  echo
  exit 1
  }
echo

# Filtering for unique transcripts
echo
echo "Filtering for unique transcripts:"
{ cut -f 10 $BLATOUT | uniq -c | awk '{print$1}' | sort | uniq -c; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Filtering of BLAT output failed. Aborting."
  echo "Check if files $INPUTFILE and $BLATOUT are correct."
  echo
  exit 1
  }
echo
echo "Filtered transcripts saved for possible later usage as"
echo "  $BLATOUT for possible later usage."
confirmgo

# Make a list of these unique transcripts (names and sequences) and convert this file to FASTA
echo
echo "Making list of unique transcripts"
cut -f 10 $BLATOUT | uniq -c | awk '{if($1==1){print $0}}' | awk '{print $2}' | awk '{printf "%012d\n", $0;}' > $UNIQUELIST || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Making list of unique transcripts failed. Aborting."
  echo "Check if files $BLATOUT and $INPUTFILE are correct."
  echo
  exit 1
  }
echo

# In order to use the join command the original transcriptome file has to be converted to TXT, the transcript numbers have to be adjusted and the file sorted
echo
echo "Converting original data into TXT for subsequent joining"
fasta2tab $INPUTFILE $INPUTTAB || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of $INPUTFILE failed. Aborting."
  echo "Check if $INPUTFILE is valid FASTA file."
  echo
  exit 1
  }
echo

echo "Sorting unique transcripts"
{ awk '{$1=sprintf("%012d", $1); print $0}' $INPUTTAB | sort > $SORTEDINPUT; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Sorting of unique transcripts failed. Aborting. Check if"
  echo "$INPUTFILE is correct FASTA file and check if file $INPUTTAB is correct."
  echo
  exit 1
  }
echo

# Apply the join command
join -j1 $UNIQUELIST $SORTEDINPUT > $JOINEDTS || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Joining failed. Aborting. Check files"
  echo "$UNIQUELIST and $SORTEDINPUT if they have same number of lines."
  echo
  exit 1
  }

# Convert to FASTA
echo
echo "Converting to FASTA"
sed 's/ /\t/g' $JOINEDTS > $JOINEDTABS || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of $JOINEDTS to FASTA failed."
  echo "Aborting. Check if file $JOINEDTS is correct."
  echo
  exit 1
  }

echo
awk '{print ">"$1"\n"$2}' $JOINEDTABS > $JOINEDFA || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of $JOINEDTS to FASTA failed."
  echo "Aborting. Check if file $JOINEDTABS is correct."
  echo
  exit 1
  }   
echo "Joined transcripts written in FASTA format as"
echo "  $JOINEDFA for possible later usage."
confirmgo

# Step 2: Find genome skimming data (only nuclear reads) which align to the unique transcripts

echo
echo "Step 2 of the pipeline - removal of reads of plastid origin."

# Get rid of the chloroplast and mitochondrial reads in the genome skimming data

# Chloroplast reads

# Create a reference plastome index with bowtie2-build
echo
echo "Creating a reference plastome index"
echo
bowtie2-build $REFERENCECP $REFERENCECP2 || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Creating a reference plastome index with bowtie2-build failed. Aborting."
  echo "Check if file $REFERENCECP is correct."
  echo
  exit 1
  }
echo

# Map the cpDNA reads to the reference plastome with BOWTIE2
echo
echo "Mapping cpDNA reads to the reference plastome. This may take longer time."
echo
bowtie2 -x $REFERENCECP2 -1 $INPUTFQ1 -2 $INPUTFQ2 -S $BOWTIE2CP || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Mapping cpDNA reads to the reference plastome with bowtie2 failed."
  echo "Aborting. Check if files $REFERENCECP,"
  echo "$REFERENCECP2, $INPUTFQ1"
  echo "and $INPUTFQ2 are correct."
  echo
  exit 1
  }
echo
echo "Mapping finished"

# Convert SAM to BAM with SAMtools
echo
echo "Converting SAM to BAM. This may take longer time."
samtools view -bT $REFERENCECP $BOWTIE2CP > $CPBAM || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of SAM to BAM with samtools failed. Aborting."
  echo "Check if files $REFERENCECP and $BOWTIE2CP are correct."
  echo
  exit 1
  }
echo "Converted file saved as"
echo "  $CPBAM for possible later usage."
confirmgo

# Remove the cpDNA reads with bam2fastq
echo
echo "Removing cpDNA reads. This may take longer time."
bam2fastq --no-aligned $CPBAM -o $FASTQNOCP#.fq || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Removal of cpDNA reads with bam2fastq failed. Aborting."
  echo "Check if file $CPBAM is correct."
  echo
  exit 1
  }
echo
echo "Removed reads saved for possible later usage as"
ls -1 $FASTQNOCP*
confirmgo

# Mitochondrial reads - optional step

if [ -n "$REFERENCEMT" ]; then

  echo
  echo "Step 3 of the pipeline - removal of reads of mitochondrial origin (optional)."

  # Create a reference mitochondrion index with bowtie2-build
  echo
  echo "Creating a reference mitochondrion index"
  echo
  bowtie2-build $REFERENCEMT $REFERENCEMT2 || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} Creating of reference mitochondrion index with bowtie2-build failed."
    echo "Aborting. Check if files $REFERENCEMT and $REFERENCEMT2 are correct."
    echo
    exit 1
    }
  echo

  # Map the mtDNA reads to the reference mitochondrion with BOWTIE2
  echo
  echo "Mapping mtDNA reads to reference mitochondrion. This may take longer time."
  echo
  bowtie2 -x $REFERENCEMT2 -1 `echo $FASTQNOCP`_1.fq -2 `echo $FASTQNOCP`_2.fq -S $BOWTIE2MT || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} Mapping mtDNA reads to reference mitochondrion with bowtie2 failed."
    echo "Aborting. Check if files $REFERENCEMT2,"
    echo "`echo $FASTQNOCP`_1.fq and `echo $FASTQNOCP`_2.fq are correct."
    echo
    exit 1
    }
  echo

  # Convert SAM to BAM with SAMtools
  echo
  echo "Converting SAM to BAM. This may take longer time."
  samtools view -bT $REFERENCEMT $BOWTIE2MT > $MTBAM || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} Conversion of SAM to BAM failed. Aborting. Check if files"
    echo "$REFERENCEMT and $BOWTIE2MT are correct."
    echo
    exit 1
    }
  echo "Converted file saved as"
  echo "  $MTBAM for possible later usage."
  confirmgo

  # Remove the mtDNA reads with bam2fastq
  echo
  echo "Removing mtDNA reads. This may take longer time."
  bam2fastq --no-aligned $MTBAM -o $FASTQNOMT#.fq || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} Removal of mtDNA reads failed. Aborting."
    echo "Check if file $MTBAM is correct."
    echo
    exit 1
    }
  echo

  # Combine the paired-end reads with FLASH - with Mt reads
  echo
  echo "Step 4 of the pipeline - combination of paired-end reads."
  echo
  echo "Combining paired-end reads"
  echo
  flash -o $FLASHOUT -M $FLASHM `echo $FASTQNOMT`_1.fq `echo $FASTQNOMT`_2.fq || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} Combining paired-end reads failed. Aborting. Check if files"
    echo "$REFERENCEMT, `echo $FASTQNOMT`_1.fq and `echo $FASTQNOMT`_2.fq are correct."
    echo
    exit 1
    }
  echo
  else
    # Combine the paired-end reads with FLASH - without Mt reads
    echo
    echo "Step 4 of the pipeline - combination of paired-end reads."
    echo
    echo "Combining paired-end reads"
    echo
    flash -o $FLASHOUT -M $FLASHM `echo $FASTQNOCP`_1.fq `echo $FASTQNOCP`_2.fq || {
      echo
      echo "${REDF}${BOLD}Error!${NORM} Combining paired-end reads failed. Aborting. Check if files"
      echo "$REFERENCECP, `echo $FASTQNOCP`_1.fq and `echo $FASTQNOCP`_2.fq are correct."
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
fastq_to_fasta -i $FLASHOUT.extendedFrags.fastq -o $FLASHOUT.extendedFrags.fa || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of FASTQ to FASTA failed. Aborting. Check if file"
  echo "$FLASHOUT.extendedFrags.fastq is correct."
  echo
  exit 1
  }
echo "Converted FASTA saved as"
echo "  $FLASHOUT.extendedFrags.fa for possible later usage."
confirmgo

# BLAT between the unique transcripts and the genome skimming data
echo
echo "BLAT between the unique transcripts and the genome skimming data."
echo "This may take longer time."
blat -t=dna -q=dna -minIdentity=$BLATIDENT -out=pslx $JOINEDFA $FLASHOUT.extendedFrags.fa $BLATOUTFIN || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} BLAT between the unique transcripts and the genome skimming data failed."
  echo "Aborting. Check if files $JOINEDFA and $FLASHOUT.extendedFrags.fa are correct."
  echo
  exit 1
  }
echo "Output saved as"
echo "  $BLATOUTFIN for possible later usage."
confirmgo

# Step 6: Assemble the obtained sequences in contigs

echo
echo "Step 6 of the pipeline - filtering of BLAT output."

# Modification of the PSLX file is needed (remove headers, select the field with the transcript (target) sequence names and the field with the query sequences, convert to FASTA) for usage in Geneious
echo
echo "Modifying PSLX BLAT output for usage in Geneious"
{ sed 1,5d $BLATOUTFIN | cut -f14,$PSLXCUT | awk '{n=split($2,a,",");for(i=1;i<=n;i++)print $1"_"NR"_"i,a[i]}' | sed 's/^/>/' | sed 's/ /\n/' > $BLATOUTFIN2; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Modifying PSLX BLAT output failed. Aborting."
  echo "Check if file $BLATOUTFIN is correct."
  echo
  exit 1
  }
echo
echo "Modified file saved as"
echo "  $BLATOUTFIN2 for possible later usage."
confirmgo

# Convert FASTA to TAB
echo
echo "Converting FASTA to TAB"
fasta2tab $BLATOUTFIN2 $TAB || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of FASTA to TAB failed. Aborting."
  echo "Check if file $BLATOUTFIN2 is correct."
  echo
  exit 1
  }
echo

{ awk '{print $1"\t"length($2)"\t"$2}' $TAB | awk '{sum+=$2}END{print sum}'; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Conversion of FASTA to TAB failed. Aborting."
  echo "Check if file $TAB is correct."
  echo
  exit 1
  }

# Remove transcripts with >1000 BLAT scores (or another value selected by user; very likely that these are repetitive elements)

# Count the number of times each transcript hit a genome skim read
echo
echo "Counting number of times each transcript hit a genom skim read"
{ cut -f1 -d_ $TAB | sort | uniq -c | sort -n -r > $TABLIST; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Counting of number of times each transcript hit a genom skim read failed."
  echo "Aborting. Check if file $TAB is correct."
  echo
  exit 1
  }

# List of the transcripts with >1000 BLAT scores (or another value selected by user)
echo
echo "Listing transcripts with >$BLATSCORE BLAT scores"
{ awk '$1>'"$BLATSCORE"'' $TABLIST | awk '{print $2}' > $TABBLAT; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Listing of transcripts with >$BLATSCORE BLAT scores failed. Aborting."
  echo "Check if file $TABLIST is correct."
  echo
  exit 1
  }

# Make a new TAB file without these transcripts
echo
echo "Removing transcripts with >$BLATSCORE BLAT score"
grep -v -f $TABBLAT $TAB > $TABREMOVED || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
  echo "Check if files $TABBLAT and $TAB are correct."
  echo
  exit 1
  }
echo

{ awk '{print $1"\t"length($2)"\t"$2}' $TABREMOVED | awk '{sum+=$2}END{print sum}'; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
  echo
  exit 1
  }

# Convert this TAB file to FASTA and remove the sequences with 'n'
echo
echo "Converting TAB to FASTA and removing sequences with \"n\""
{ grep -v n $TABREMOVED | sed 's/^/>/' | sed 's/\t/\n/' > $FINALA; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
  echo "Check if file $TABREMOVED is correct."
  echo
  exit 1
  }

grep -v n $TABREMOVED | awk '{print $1"\t"length($2)}' | awk '{s+=$2;a++}END{print s}' || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting."
  echo "Check if file $TABREMOVED is correct."
  echo
  exit 1
  }

# Remove unneeded temporal files - keep only *.pslx, *.fasta and *.bam
echo "Removing unneeded temporal files"
if [ -n "$REFERENCEMT" ]; then
  rm $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $REFERENCEMT2* $BOWTIE2MT $FLASHOUT.extendedFrags.fastq $TAB $TABLIST $TABBLAT $TABREMOVED || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} Removal of temporal files failed. Remove following files manually:"
    echo "$UNIQUELIST, $INPUTTAB, $SORTEDINPUT,"
    echo "$JOINEDTS, $JOINEDTABS, $REFERENCECP2*,"
    echo "$BOWTIE2CP, $REFERENCEMT2*, $BOWTIE2MT,"
    echo "$FLASHOUT.extendedFrags.fastq, $TAB, $TABLIST,"
    echo "$TABBLAT and $TABREMOVED."
    confirmgo
    }
  else
    rm $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $FLASHOUT.extendedFrags.fastq $TAB $TABLIST $TABBLAT $TABREMOVED || {
      echo
      echo "${REDF}${BOLD}Error!${NORM} Removal of temporal files failed. Remove following files manually:"
      echo "$UNIQUELIST, $INPUTTAB, $SORTEDINPUT,"
      echo "$JOINEDTS, $JOINEDTABS, $REFERENCECP2*,"
      echo "$BOWTIE2CP, $FLASHOUT.extendedFrags.fastq, $TAB,"
      echo "$TABLIST, $TABBLAT and $TABREMOVED."
      confirmgo
      }
  fi

# List kept files which user can use for another analysis
echo
echo "Following files are kept for possible later usage (see manual for details):"
echo "${BLUF}================================================================================${NORM}"
if [ -n "$REFERENCEMT" ]; then
  if [ "$RENAMEDFASTASEQENCES" == "YES" ]; then
    echo "0)  List of old names of FASTA sequences in $INPUTFILE0 and"
    echo "    renamed FASTA sequences in $INPUTFILE:"
    echo "$TRANSCRIPTOMEFASTANAMES"
    fi
  echo "1)  Output of BLAT (removal of transcripts sharing ≥90% sequence similarity):"
  echo "$BLATOUT"
  echo "2)  Unique transcripts in FASTA format:"
  echo "$JOINEDFA"
  echo "3)  SAM converted to BAM (removal of reads of plastid origin):"
  echo "$CPBAM"
  echo "4)  Genome skim data without cpDNA reads:"
  ls $FASTQNOCP*
  echo "5)  SAM converted to BAM (removal of reads of mitochondrial origin):"
  echo "$MTBAM"
  echo "6)  Genome skim data without mtDNA reads:"
  ls $FASTQNOMT*
  echo "7)  Combined paired-end genome skim reads:"
  echo "$FLASHOUT.extendedFrags.fa"
  echo "8)  Output of BLAT (matching of the unique transcripts and the filtered,"
  echo "    combined genome skim reads sharing ≥85% sequence similarity):"
  echo "$BLATOUTFIN"
  echo "9)  Matching sequences in FASTA:"
  echo "$BLATOUTFIN2"
  echo "10) Final FASTA sequences for usage in Geneious:"
  echo "$FINALA"
  else
    if [ "$RENAMEDFASTASEQENCES" == "YES" ]; then
      echo "0)  List of old names of FASTA sequences in $INPUTFILE0 and"
      echo "    renamed FASTA sequences in $INPUTFILE:"
      echo "$TRANSCRIPTOMEFASTANAMES"
      fi
    echo "1) Output of BLAT (removal of transcripts sharing ≥90% sequence similarity):"
    echo "$BLATOUT"
    echo "2) Unique transcripts in FASTA format:"
    echo "$JOINEDFA"
    echo "3) SAM converted to BAM (removal of reads of plastid origin):"
    echo "$CPBAM"
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
echo "${BLUF}================================================================================${NORM}"
confirmgo

echo
echo "Success!"
echo
echo "Resulting FASTA was saved as"
echo "${BOLD}$FINALA${NORM} for usage in Geneious (step 7 of the pipeline)."
echo "Use this file in next step of the pipeline. See README and manual for details."
echo "${BLUF}================================================================================${NORM}"
confirmgo

echo "Run Geneious (tested with versions 6, 7 and 8). See README and manual for details."
echo
echo "Import output file \"$FINALA\" (File | Import | From File...)."
echo "Select the file and go to menu Tools | Align / Assemble | De Novo Assemble."
echo "In \"Data\" frame select \"Assembe by 1st (...) Underscore\"."
echo "In \"Method\" frame select Geneious Assembler (if you don't have other assemblers,"	
echo "  this option might be missing) and \"Medium Sensitivity / Fast\" Sensitivity."
echo "In \"Results\" frame check \"Save assembly report\", \"Save list of unused reads\","
echo "  \"Save in sub-folder\", \"Save contigs\" (do not check \"Maximum\") and"
echo "  \"Save consensus sequences\"."
echo "Do not trim. Otherwise keep defaults. Run it."
echo "Geneious may warn about possible hanging because of big file."
echo "Do not use Geneious for other tasks during the assembly."
echo
echo "Select all resulting contigs (typically named \"* Contig #\") and export them"
echo "  (File | Export | Selected Documents...) as \"Tab-separated table values (*.tsv)\"."
echo "Save following columns (Hold Ctrl key to mark more fields): \"# Sequences\","
echo "  \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and"
echo "  \"Sequence Length\"."
echo "If this option would be inaccessible for you, export all columns."
echo "Warning! Do not select and export \"* Consensus Sequences\", \"* Unused Reads\""
echo "or \"* Assembly Report\" - only the individual  \"* contig #\" files."
echo
echo "Select items \"Consensus Sequences\" and \"Unused Reads\" and export both of them"
echo "  as one FASTA file."
echo "Go to menu File | Export | Selected Documents... and choose FASTA file type."
echo
echo "Use exported files from Geneious as input for part B of the script"
echo "  (sondovac_part_b.sh; see README and/or PDF manual for details)."
echo
echo "Script exited successfully..."
echo

exit
