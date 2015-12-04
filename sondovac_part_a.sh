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
      echo -e "\t${REDF}-a${NORM}\t${CYAF}Maximum overlap length expected in approximately 90% of read"
      echo -e "\t\t  pairs${NORM} (parameter \"-M\" of FLASH, see its manual for details)."
      echo -e "\t\tDefault value: 65 (integer ranging from 10 to 300)"
      echo -e "\t${REDF}-y${NORM}\t${CYAF}Sequence similarity between unique transcripts and the filtered,"
      echo -e "\t\t  combined genome skim reads${NORM} (parameter \"-minIdentity\" of BLAT,"
      echo -e "\t\t  see its manual for details)."
      echo -e "\t\tDefault value: 85 (integer ranging from 70 to 100; the default"
      echo -e "\t\t  value of 85% minimum sequence similarity suggests gene"
      echo -e "\t\t  orthology)"
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
      echo "${CYAF}Running in interactive mode..${NORM}."
      STARTINI="I"
      CHECKMODE=$((CHECKMODE+1))
      ;;
    n)
      echo "${CYAF}Running in non-interactive mode...${NORM}"
      STARTINI="N"
      CHECKMODE=$((CHECKMODE+1))
      ;;
    f)
      INPUTFILE0=$OPTARG
      echo "Transcriptome file: ${REDF}$INPUTFILE0${NORM}"
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
      if [[ "$FLASM" =~ ^[0-9]+$ ]] && [ "$FLASM" -ge 10 -a "$FLASM" -le 300 ]; then
	echo "Maximum overlap length expected in approximately 90% of read pairs: ${REDF}$FLASHM${NORM}"
	else
	  echo "${REDF}${BOLD}Error!${NORM} For parameter \"-a\" you did not provide an integer ranging from 10 to 300!"
	  echo
	  exit 1
	fi
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

# Set variables for working directory and ${BOLD}PATH${NORM}
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
    echo "  download bowtie2-$BOWTIE2V-source.zip, compile it and ensure it is in ${BOLD}PATH${NORM}."
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
  echo "  and \"bowtie2-build-s\" are required but not installed or available in ${BOLD}PATH${NORM}.${NORM}"
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile \"Bowtie2-$BOWTIE2V\" from source available together with this script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to compile \"Bowtie2-$BOWTIE2V\" from source code${NORM} downloaded from"
    echo "  ${REDF}http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE2V/${NORM}"
    echo "Type \"${REDF}D${NORM}\" ${CYAF}to download \"Bowtie2-$BOWTIE2V\" binary${NORM} from"
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
	    echo "Downloading \"${REDF}Bowtie2-$BOWTIE2V${NORM}\" binaries for $OS $OSB" &&
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
	    echo "Downloading \"${REDF}Bowtie2-$BOWTIE2V${NORM}\" binaries for $OS $OSB" &&
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
	    echo "Downloading \"${REDF}Bowtie2-$BOWTIE2V${NORM}\" binaries for $OS $OSB" &&
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
	    echo "Unknown OS or OS without \"${REDF}Bowtie2-$BOWTIE2V${NORM}\" binary available."
	    compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
	  fi
	  break
	  ;;
	S|s)
	  echo
	  downloaderselector
	  checktools unzip
	  echo
	  echo "Downloading \"${REDF}Bowtie2-$BOWTIE2V${NORM}\" source code"
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
	      echo "Binary is not available for ${REDF}$OS $OSB${NORM}. Going to compile it from source code."
	      echo
	      compilebowtie $SCRIPTDIR/src/bowtie2-$BOWTIE2V
	      ;;
	  esac
	  break
	  ;;
	H|h)
	  if [ "$OS" == "Mac" ]; then
	    { echo "Installing \"${REDF}Bowtie2-$BOWTIE2V${NORM}\" using Homebrew" &&
	    brew install homebrew/science/bowtie2 &&
	    echo "\"${REDF}Bowtie2${NORM}\" is available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of \"Bowtie2\" failed.${NORM} Please, do it manually. For details see"
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
	  echo " latest Bowtie2 and ensure it is in ${BOLD}PATH${NORM}."
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
  echo >&2 "\"${REDF}samtools${NORM}\" is required but not installed or available in ${BOLD}PATH${NORM}."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile \"SAMtools-1.2\" from source available together with this script.${NORM}"
    echo "  Makefile was modified not to require GNU ncurses library."
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download latest developmental \"SAMtools\" source${NORM} from"
    echo "  ${REDF}https://github.com/samtools/samtools/${NORM} and compile it. Compilation requires GNU"
    echo "  ncurses library and is recommended only for advanced users. If compilation"
    echo "  fails, check SAMtools' ${REDF}INSTALL${NORM} file for details and adjust its Makefile."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy \"SAMtools-1.2\" binary${NORM} available together with the script"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
    echo "  See \"${REDF}brew info homebrew/science/samtools${NORM}\" for more details."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install \"SAMtools\" yourselves."
    read SAMTOOLS
    while :
    do
      case "$SAMTOOLS" in
	C|c)
	  compilesamtools || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to ${REDF}http://www.htslib.org/download/${NORM}"
	    echo "  download samtools-1.2, compile it and ensure it is in ${BOLD}PATH${NORM}."
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
	  echo "Downloading \"${REDF}SAMtools${NORM}\" sources..." &&
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
	    echo "  latest stable version of \"${REDF}SAMtools${NORM}\" or ${REDF}https://github.com/samtools/${NORM} for latest"
	    echo "  developmental version, download samtools, compile it and ensure it is in ${BOLD}PATH${NORM}."
	    echo
	    exit 1
	    }
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying \"${REDF}SAMtools${NORM}\" binaries"
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
	      echo "Copying \"${REDF}SAMtools${NORM}\" binaries"
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
	      echo "Binary is not available for ${REDF}$OS $OSB${NORM}."
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
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of \"SAMtools\" failed.${NORM} Please, do it manually. For details see"
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
	  echo "Please, go to ${REDF}http://www.htslib.org/${NORM} and install \"${REDF}SAMtools${NORM}\" and ensure it is"
	  echo "  in ${BOLD}PATH${NORM}."
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
  echo "Compiling \"${REDF}bam2fastq${NORM}\"..." &&
  cd $1 &&
  make -s &&
  cp bam2fastq $BIN/ &&
  cd $WORKDIR &&
  echo "\"${REDF}bam2fastq${NORM}\" is available. ${GREF}OK.${NORM}"
  }

# Check if bam2fastq is available
{ command -v bam2fastq >/dev/null 2>&1 && echo "\"${REDF}bam2fastq${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "\"${REDF}bam2fastq${NORM}\" is required but not installed or available in ${BOLD}PATH${NORM}."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile \"bam2fastq\" from source available together with this script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download \"bam2fastq\" source${NORM} from"
    echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM} and compile it."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy \"bam2fastq\" binary available together with the script${NORM}"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install \"bam2fastq\" yourselves."
    read bam2fastq
    while :
    do
      case "$bam2fastq" in
	C|c)
	  compilebam2fastq $SCRIPTDIR/src/bam2fastq-1.1.0 || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
	    echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}, download bam2fastq,"
	    echo "  compile it manually and ensure it is in ${BOLD}PATH${NORM}."
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
	    echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}, download \"${REDF}bam2fastq${NORM}\","
	    echo "  compile it manually and ensure it is in ${BOLD}PATH${NORM}."
	    echo
	    exit 1
	    }
		break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying \"${REDF}bam2fastq${NORM}\" binary"
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/bam2fastq $BIN/
	      break
	      ;;
	    Linux)
	      echo "Copying \"${REDF}bam2fastq${NORM}\" binary"
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/bam2fastq $BIN/
	      break
	      ;;
	    *) echo
	      echo "Binary is not available for ${REDF}$OS $OSB${NORM}."
	      echo
	      compilebam2fastq $SCRIPTDIR/src/bam2fastq-1.1.0 || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Compilation failed.${NORM} Please, go to"
		echo "  ${REDF}http://gsl.hudsonalpha.org/information/software/bam2fastq${NORM}, download \"${REDF}bam2fastq${NORM}\","
		echo "  compile it manually and ensure it is in ${BOLD}PATH${NORM}."
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
	  echo "  download \"${REDF}bam2fastq${NORM}\", compile it manually and ensure it is in ${BOLD}PATH${NORM}."
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
  echo "Compiling \"${REDF}FLASH${NORM}\" from source code..." &&
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
    echo "  FLASH-*.tar.gz, compile it manually and ensure it is in ${BOLD}PATH${NORM}."
    echo
    exit 1
    }
  }

# Check if FLASH is available
{ command -v flash >/dev/null 2>&1 && echo "\"${REDF}flash${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "${CYAF}FLASH is required but not installed or available in ${BOLD}PATH${NORM}.${NORM}"
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile \"FLASH\" from source available together with this script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download \"FLASH\" source${NORM} from"
    echo "  ${REDF}http://sourceforge.net/projects/flashpage/${NORM} and compile it."
    echo "Type \"${REDF}D${NORM}\" ${CYAF}to download \"FLASH\" binary${NORM} from"
    echo "  ${REDF}http://sourceforge.net/projects/flashpage/${NORM} (available only for Windows)."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy \"FLASH 1.2.11\" binary available together with the script${NORM}"
    echo "  (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
    echo "  See \"${REDF}brew info homebrew/science/flash${NORM}\" for more details."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install \"FLASH\" yourselves."
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
	  echo "Downloading \"${REDF}FLASH${NORM}\" source code"
	  $DOWNLOADER FLASH-1.2.11.tar.gz http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz || {
	    echo
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to ${REDF}http://sourceforge.net/projects/flashpage/files/${NORM}"
	    echo "  download FLASH-*.tar.gz, compile it manually and ensure it is in ${BOLD}PATH${NORM}."
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
	      echo "Downloading \"${REDF}FLASH${NORM}\" for $OS"
	      $DOWNLOADER FLASH-1.2.11-windows-bin.zip http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11-windows-bin.zip || {
		echo
		echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to ${REDF}http://sourceforge.net/projects/flashpage/files/${NORM}"
		echo "  download FLASH-*windows-bin.zip, unpack it and ensure it is in ${BOLD}PATH${NORM}."
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
	      echo "Copying \"${REDF}FLASH${NORM}\" binary"
	      cp -p $SCRIPTDIR/pkgs/macosx/bin/flash $BIN/
	      break
	      ;;
	    Linux)
	      echo "Copying \"${REDF}FLASH${NORM}\" binary"
	      cp -p $SCRIPTDIR/pkgs/linux64b/bin/flash $BIN/
	      break
	      ;;
	    *) echo
	      echo "Binary is not available for ${REDF}$OS $OSB${NORM}."
	      echo
	      compileflash $SCRIPTDIR/src/FLASH-1.2.11
	      ;;
	  esac
	  break
	  ;;
	H|h)
	  if [ "$OS" == "Mac" ]; then			
	    { echo "Installing \"${REDF}FLASH${NORM}\" using Homebrew" &&
	    brew install homebrew/science/flash &&
	    echo "\"${REDF}FLASH${NORM}\" is available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of \"FLASH\" failed.${NORM} Please, do it manually. For details see"
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
	  echo "Please, go to ${REDF}http://ccb.jhu.edu/software/FLASH/${NORM} and install \"${REDF}FLASH${NORM}\" and ensure it"
	  echo "  is in ${BOLD}PATH${NORM}"
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
  echo "Compiling \"${REDF}FASTX-Toolkit${NORM}\" from source code..." &&
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
    echo "  they are in ${BOLD}PATH${NORM}."
    echo
    exit 1
    }
  }

# Check if fastq_to_fasta is available
{ command -v fastq_to_fasta >/dev/null 2>&1 && echo "\"${REDF}fastq_to_fasta${NORM}\" is available. ${GREF}OK.${NORM}"; } || {
  echo
  echo >&2 "\"${REDF}fastq_to_fasta${NORM}\" is required but not installed or available in ${BOLD}PATH${NORM}."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"${REDF}C${NORM}\" ${CYAF}to compile \"FASTX-Toolkit\" from source available together with this"
    echo "  script.${NORM}"
    echo "Type \"${REDF}S${NORM}\" ${CYAF}to download \"FASTX-Toolkit\" source${NORM} from"
    echo "  ${REDF}http://hannonlab.cshl.edu/fastx_toolkit/${NORM} and compile it."
    echo "Type \"${REDF}B${NORM}\" ${CYAF}to copy \"FASTX-Toolkit\" 0.0.14 binary${NORM} available together with the"
    echo "  script (recommended, available for Linux and Mac OS X)."
    echo "Type \"${REDF}H${NORM}\" ${CYAF}for installation using Homebrew${NORM} (only for Mac OS X, recommended)."
    echo "  See \"${REDF}brew info homebrew/science/fastx_toolkit${NORM}\" for more details."
    echo "Type \"${REDF}M${NORM}\" ${CYAF}for manual installation${NORM} - script will exit and you will have to"
    echo "  install \"FASTX-Toolkit\" yourselves."
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
	    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Download failed.${NORM} Please, go to"
	    echo "  ${REDF}https://github.com/agordon/fastx_toolkit/releases/download/${NORM} and download"
	    echo "  libgtextutils and fastx_toolkit, install them and ensure they are in ${BOLD}PATH${NORM}."
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
	      echo "Copying \"${REDF}FASTX-Toolkit${NORM}\" binaries"
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
	      echo "Copying \"${REDF}FASTX-Toolkit${NORM}\" binaries"
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
	      echo "Binary is not available for ${REDF}$OS $OSB${NORM}."
	      echo
	      compilefastx $SCRIPTDIR/src/libgtextutils-0.7/ $SCRIPTDIR/src/fastx_toolkit-0.0.14/
	      ;;
	  esac
	  break
	  ;;
	H|h)
	  if [ "$OS" == "Mac" ]; then
	    { echo "Installing \"${REDF}FASTX-Toolkit${NORM}\" using Homebrew" &&
	    brew install homebrew/science/fastx_toolkit &&
	    echo "\"${REDF}FASTX-Toolkit${NORM}\" is available. ${GREF}OK.${NORM}"
	    } || {
	      echo
	      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Installation of \"FASTX-Toolkit\" failed.${NORM} Please, do it manually. For details"
	      echo "see \"brew info ${REDF}homebrew/science/fastx_toolkit${NORM}\" and \"${REDF}brew help${NORM}\"."
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
	  echo "Please, go to ${REDF}http://hannonlab.cshl.edu/fastx_toolkit/download.html${NORM} download"
	  echo "  libgtextutils-*.tar.gz and fastx_toolkit-*.tar.bz2, compile them and ensure"
	  echo "  they are in ${BOLD}PATH${NORM}"
	  exit 2
	  ;;
	*) echo "Wrong option. Use C, S, B, H or M." && read FASTX;;
      esac
    done
  else
    exit 1
  fi
  }

echo
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
  echo "${CYAF}Would you like to use mitochondriome reference${NORM} sequence input file in FASTA"
  echo "  format? (${REDF}Yes${NORM}/${REDF}No${NORM})"
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
	echo "OK, we will not use mitochondriome reference sequence. ${CYAF}Continuing.${NORM}"
	REFERENCEMT=""
	break
	;;
      *) echo "${CYAF}Wrong option.${NORM} Use ${REDF}Y${NORM} or ${REDF}N${NORM}." && read $MTINPUTQ;;
    esac
  done
  echo
  fi

if [ -z "$REFERENCEMT" ]; then
  echo
  echo "${CYAF}${BOLD}Warning!${NORM} There is no mitochondriome reference sequence."
  echo
  else
    echo
  fi

# Input file in FASTA format
echo "Input file: ${REDF}$INPUTFILE0${NORM}"
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
echo "Input file: ${REDF}$REFERENCECP${NORM}"
# Reference genome - plastome index - temporary file - will be deleted
REFERENCECP2="${OUTPUTFILENAME%.*}.cp"
# Input cpDNA reads in FASTQ
echo "Input file: ${REDF}$INPUTFQ1${NORM}"
echo "Input file: ${REDF}$INPUTFQ2${NORM}"
# cpDNA reads mapped to reference - temporary file - will be deleted
BOWTIE2CP="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads.sam"
# SAM converted into BAM (removal of reads of plastid origin)
CPBAM="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads.bam"
# Genome skim data without cpDNA reads
FASTQNOCP="${OUTPUTFILENAME%.*}_genome_skim_data_no_cp_reads"
# Input - reference genome - mtDNA
if [ -n "$REFERENCEMT" ]; then
  echo "Input file: ${REDF}$REFERENCEMT${NORM}"
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
  echo "Input file ${REDF}$1${NORM}"
  echo "  ${CYAF}contains some non-allowed names.${NORM} New file with correct names"
  echo "  (only increasing unique numbers) will be created."
  echo "Depending on size of your transcriptome file ${CYAF}it can take longer time.${NORM}"
  while read LINE; do
    N=$((++N)) &&
    echo $LINE | sed -e 's/^>.*$/>'$N'/'
    done < $1 > $2 || {
      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Renaming failed.${NORM} Aborting. Ensure ${REDF}$1${NORM}"
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
  echo "File ${REDF}$TRANSCRIPTOMEFASTANAMES${NORM} contains two columns:"
  echo "  ${CYAF}1)${NORM} Old names in original ${REDF}$1${NORM} and"
  echo "  ${CYAF}2)${NORM} New names in required format as in ${REDF}$2${NORM}"
  echo "  The sequences remain intact. This might be needed to trace back some sequences."
  confirmgo
  RENAMEDFASTASEQENCES="YES"
  }

echo
echo "${CYAF}Checking if transcriptome input file ${REDF}$INPUTFILE0${NORM}"
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
      echo "Input file ${REDF}$INPUTFILE0${NORM}"
      echo " has correct names. ${CYAF}Continuing.${NORM}"
      INPUTFILE=$INPUTFILE0
      RENAMEDFASTASEQENCES="NO"
    fi

# Step 1: Obtain unique transcripts.

echo
echo "${REDF}Step 1 of the pipeline${NORM} - removal of transcripts sharing ≥90% sequence similarity."

# BLAT between the transcriptome itself
echo
echo "${CYAF}Launching BLAT between the transcriptome itself${NORM}"
blat -t=dna -q=dna -noHead -out=psl $INPUTFILE $INPUTFILE $BLATOUT || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}BLAT failed.${NORM} Aborting. Check if ${REDF}$INPUTFILE${NORM} is correct."
  echo
  exit 1
  }

# Filtering for unique transcripts
echo
echo "${CYAF}Filtering for unique transcripts:${NORM}"
{ cut -f 10 $BLATOUT | uniq -c | awk '{print$1}' | sort | uniq -c; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Filtering of BLAT output failed.${NORM} Aborting."
  echo "Check if files ${REDF}$INPUTFILE${NORM} and ${REDF}$BLATOUT${NORM} are correct."
  echo
  exit 1
  }
echo
echo "${CYAF}Filtered transcripts saved${NORM}  for possible later usage as"
echo "  ${REDF}$BLATOUT${NORM}  for possible later usage."
confirmgo

# Make a list of these unique transcripts (names and sequences) and convert this file to FASTA
echo "Making list of unique transcripts"
cut -f 10 $BLATOUT | uniq -c | awk '{if($1==1){print $0}}' | awk '{print $2}' | awk '{printf "%012d\n", $0;}' > $UNIQUELIST || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Making list of unique transcripts failed.${NORM} Aborting."
  echo "Check if files ${REDF}$BLATOUT${NORM} and ${REDF}$INPUTFILE${NORM} are correct."
  echo
  exit 1
  }

# In order to use the join command the original transcriptome file has to be converted to TXT, the transcript numbers have to be adjusted and the file sorted
echo
echo "Converting original data into TXT for subsequent joining"
fasta2tab $INPUTFILE $INPUTTAB || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of${NORM} ${REDF}$INPUTFILE${NORM} ${CYAF}failed.${NORM} Aborting."
  echo "Check if ${REDF}$INPUTFILE${NORM} is valid FASTA file."
  echo
  exit 1
  }
echo

echo "${CYAF}Sorting unique transcripts${NORM}"
{ awk '{$1=sprintf("%012d", $1); print $0}' $INPUTTAB | sort > $SORTEDINPUT; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Sorting of unique transcripts failed.${NORM} Aborting. Check if"
  echo "${REDF}$INPUTFILE${NORM} is correct FASTA file and check if file ${REDF}$INPUTTAB${NORM} is correct."
  echo
  exit 1
  }

# Apply the join command
join -j1 $UNIQUELIST $SORTEDINPUT > $JOINEDTS || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Joining failed.${NORM} Aborting. Check files"
  echo "${REDF}$UNIQUELIST${NORM} and ${REDF}$SORTEDINPUT${NORM} if they have same number of lines."
  echo
  exit 1
  }

# Convert to FASTA
echo
echo "Converting to FASTA"
sed 's/ /\t/g' $JOINEDTS > $JOINEDTABS || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of${NORM} ${REDF}$JOINEDTS${NORM} ${CYAF}to FASTA failed.${NORM}"
  echo "Aborting. Check if file  ${REDF}$JOINEDTS${NORM} is correct."
  echo
  exit 1
  }

echo
awk '{print ">"$1"\n"$2}' $JOINEDTABS > $JOINEDFA || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of${NORM} ${REDF}$JOINEDTS${NORM} ${CYAF}to FASTA failed.${NORM}"
  echo "Aborting. Check if file ${REDF}$JOINEDTABS${NORM} is correct."
  echo
  exit 1
  }
echo "${CYAF}Joined transcripts written in FASTA${NORM} format as"
echo "  ${REDF}$JOINEDFA${NORM} for possible later usage."
confirmgo

# Step 2: Find genome skimming data (only nuclear reads) which align to the unique transcripts

echo "${REDF}Step 2 of the pipeline${NORM} - removal of reads of plastid origin."

# Get rid of the chloroplast and mitochondrial reads in the genome skimming data

# Chloroplast reads

# Create a reference plastome index with bowtie2-build
echo
echo "${CYAF}Creating a reference plastome index${NORM}"
echo
bowtie2-build $REFERENCECP $REFERENCECP2 || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Creating a reference plastome index with bowtie2-build failed.${NORM}  Aborting."
  echo "Check if file ${REDF}$REFERENCECP${NORM} is correct."
  echo
  exit 1
  }
echo

# Map the cpDNA reads to the reference plastome with BOWTIE2
echo "${CYAF}Mapping cpDNA reads to the reference plastome. This may take longer time.${NORM}"
echo
bowtie2 -x $REFERENCECP2 -1 $INPUTFQ1 -2 $INPUTFQ2 -S $BOWTIE2CP || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Mapping cpDNA reads to the reference plastome with bowtie2 failed.${NORM}"
  echo "Aborting. Check if files ${REDF}$REFERENCECP${NORM},"
  echo "${REDF}$REFERENCECP2${NORM}, ${REDF}$INPUTFQ1${NORM}"
  echo "and ${REDF}$INPUTFQ2${NORM} are correct."
  echo
  exit 1
  }
echo
echo "Mapping finished"

# Convert SAM to BAM with SAMtools
echo
echo "${CYAF}Converting SAM to BAM. This may take longer time.${NORM}"
samtools view -bT $REFERENCECP $BOWTIE2CP > $CPBAM || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of SAM to BAM with samtools failed.${NORM} Aborting."
  echo "Check if files ${REDF}$REFERENCECP${NORM} and ${REDF}$BOWTIE2CP${NORM} are correct."
  echo
  exit 1
  }
echo "${CYAF}Converted BAM file saved${NORM} as"
echo "  ${REDF}$CPBAM${NORM} for possible later usage."
confirmgo

# Remove the cpDNA reads with bam2fastq
echo "${CYAF}Removing cpDNA reads. This may take longer time.${NORM}"
bam2fastq --no-aligned $CPBAM -o $FASTQNOCP#.fq || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removal of cpDNA reads with bam2fastq failed.${NORM} Aborting."
  echo "Check if file ${REDF}$CPBAM${NORM} is correct."
  echo
  exit 1
  }
echo
echo "${CYAF}Removed reads saved${NORM} for possible later usage as"
ls -1 $FASTQNOCP*
confirmgo

# Mitochondrial reads - optional step

if [ -n "$REFERENCEMT" ]; then

  echo "${REDF}Step 3 of the pipeline${NORM} - removal of reads of mitochondrial origin (optional)."

  # Create a reference mitochondrion index with bowtie2-build
  echo
  echo "${CYAF}Creating a reference mitochondrion index${NORM}"
  echo
  bowtie2-build $REFERENCEMT $REFERENCEMT2 || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Creating of reference mitochondrion index with bowtie2-build failed.${NORM}"
    echo "Aborting. Check if files ${REDF}$REFERENCEMT${NORM} and ${REDF}$REFERENCEMT2${NORM} are correct."
    echo
    exit 1
    }
  echo

  # Map the mtDNA reads to the reference mitochondrion with BOWTIE2
  echo "${CYAF}Mapping mtDNA reads to reference mitochondrion. This may take longer time.${NORM}"
  echo
  bowtie2 -x $REFERENCEMT2 -1 `echo $FASTQNOCP`_1.fq -2 `echo $FASTQNOCP`_2.fq -S $BOWTIE2MT || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Mapping mtDNA reads to reference mitochondrion with bowtie2 failed.${NORM}"
    echo "Aborting. Check if files ${REDF}$REFERENCEMT2${NORM},"
    echo "${REDF}`echo $FASTQNOCP`_1.fq${NORM} and ${REDF}`echo $FASTQNOCP`_2.fq${NORM} are correct."
    echo
    exit 1
    }
  echo

  # Convert SAM to BAM with SAMtools
  echo "${CYAF}Converting SAM to BAM. This may take longer time.${NORM}"
  samtools view -bT $REFERENCEMT $BOWTIE2MT > $MTBAM || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of SAM to BAM failed.${NORM} Aborting. Check if files"
    echo "${REDF}$REFERENCEMT${NORM} and ${REDF}$BOWTIE2MT${NORM} are correct."
    echo
    exit 1
    }
  echo "${CYAF}Converted file saved${NORM} as"
  echo "  ${REDF}$MTBAM${NORM} for possible later usage."
  confirmgo

  # Remove the mtDNA reads with bam2fastq
  echo "${CYAF}Removing mtDNA reads. This may take longer time.${NORM}"
  bam2fastq --no-aligned $MTBAM -o $FASTQNOMT#.fq || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removal of mtDNA reads failed.${NORM} Aborting."
    echo "Check if file ${REDF}$MTBAM${NORM} is correct."
    echo
    exit 1
    }
  echo

  # Combine the paired-end reads with FLASH - with Mt reads
  echo "${REDF}Step 4 of the pipeline${NORM} - combination of paired-end reads."
  echo
  echo "${CYAF}Combining paired-end reads${NORM}"
  echo
  flash -o $FLASHOUT -M $FLASHM `echo $FASTQNOMT`_1.fq `echo $FASTQNOMT`_2.fq || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Combining paired-end reads failed.${NORM} Aborting. Check if files"
    echo "${REDF}$REFERENCEMT${NORM}, ${REDF}`echo $FASTQNOMT`_1.fq${NORM} and ${REDF}`echo $FASTQNOMT`_2.fq${NORM} are correct."
    echo
    exit 1
    }
  echo
  else
    # Combine the paired-end reads with FLASH - without Mt reads
    echo "${REDF}Step 4 of the pipeline${NORM} - combination of paired-end reads."
    echo
    echo "${CYAF}Combining paired-end reads${NORM}"
    echo
    flash -o $FLASHOUT -M $FLASHM `echo $FASTQNOCP`_1.fq `echo $FASTQNOCP`_2.fq || {
      echo
      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Combining paired-end reads failed.${NORM} Aborting. Check if files"
      echo "${REDF}$REFERENCECP${NORM}, ${REDF}`echo $FASTQNOCP`_1.fq${NORM} and ${REDF}`echo $FASTQNOCP`_2.fq${NORM} are correct."
      echo
      exit 1
      }
  fi

# Step 5 - matching of the unique transcripts and the filtered, combined genome skim reads sharing ≥85% sequence similarity

echo "${REDF}Step 5 of the pipeline${NORM} - matching of the unique transcripts and the filtered,"
echo "  combined genome skim reads sharing ≥85% sequence similarity."

# Convert FASTQ file to FASTA
echo
echo "${CYAF}Converting FASTQ to FASTA. This may take longer time.${NORM}"
fastq_to_fasta -i $FLASHOUT.extendedFrags.fastq -o $FLASHOUT.extendedFrags.fa || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of FASTQ to FASTA failed.${NORM} Aborting. Check if file"
  echo "${REDF}$FLASHOUT.extendedFrags.fastq${NORM} is correct."
  echo
  exit 1
  }
echo "${CYAF}Converted FASTA saved${NORM} as"
echo "  ${REDF}$FLASHOUT.extendedFrags.fa${NORM} for possible later usage."
confirmgo

# BLAT between the unique transcripts and the genome skimming data
echo "${CYAF}BLAT between the unique transcripts and the genome skimming data."
echo "This may take longer time.${NORM}"
blat -t=dna -q=dna -minIdentity=$BLATIDENT -out=pslx $JOINEDFA $FLASHOUT.extendedFrags.fa $BLATOUTFIN || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}BLAT between the unique transcripts and the genome skimming data failed.${NORM}"
  echo "Aborting. Check if files ${REDF}$JOINEDFA${NORM} and ${REDF}$FLASHOUT.extendedFrags.fa${NORM} are correct."
  echo
  exit 1
  }
echo "${CYAF} BLAT output saved${NORM} as"
echo "  ${REDF}$BLATOUTFIN${NORM} for possible later usage."
confirmgo

# Step 6: Assemble the obtained sequences in contigs

echo "${REDF}Step 6 of the pipeline${NORM} - filtering of BLAT output."

# Modification of the PSLX file is needed (remove headers, select the field with the transcript (target) sequence names and the field with the query sequences, convert to FASTA) for usage in Geneious
echo
echo "${CYAF}Modifying PSLX BLAT output${NORM} for usage in Geneious"
{ sed 1,5d $BLATOUTFIN | cut -f14,$PSLXCUT | awk '{n=split($2,a,",");for(i=1;i<=n;i++)print $1"_"NR"_"i,a[i]}' | sed 's/^/>/' | sed 's/ /\n/' > $BLATOUTFIN2; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Modifying PSLX BLAT output failed.${NORM} Aborting."
  echo "Check if file ${REDF}$BLATOUTFIN${NORM} is correct."
  echo
  exit 1
  }
echo
echo "${CYAF}Modified file saved${NORM} as"
echo "  ${REDF}$BLATOUTFIN2${NORM} for possible later usage."
confirmgo

# Convert FASTA to TAB
echo "Converting FASTA to TAB"
fasta2tab $BLATOUTFIN2 $TAB || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of FASTA to TAB failed.${NORM} Aborting."
  echo "Check if file ${REDF}$BLATOUTFIN2${NORM} is correct."
  echo
  exit 1
  }
echo

{ awk '{print $1"\t"length($2)"\t"$2}' $TAB | awk '{sum+=$2}END{print sum}'; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Conversion of FASTA to TAB failed.${NORM} Aborting."
  echo "Check if file ${REDF}$TAB${NORM} is correct."
  echo
  exit 1
  }

# Remove transcripts with >1000 BLAT scores (or another value selected by user; very likely that these are repetitive elements)

# Count the number of times each transcript hit a genome skim read
echo
echo "Counting number of times each transcript hit a genom skim read"
{ cut -f1 -d_ $TAB | sort | uniq -c | sort -n -r > $TABLIST; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Counting of number of times each transcript hit a genom skim read failed.${NORM}"
  echo "Aborting. Check if file ${REDF}$TAB${NORM} is correct."
  echo
  exit 1
  }

# List of the transcripts with >1000 BLAT scores (or another value selected by user)
echo
echo "Listing transcripts with >${CYAF}$BLATSCORE${NORM} BLAT scores"
{ awk '$1>'"$BLATSCORE"'' $TABLIST | awk '{print $2}' > $TABBLAT; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Listing of transcripts with >$BLATSCORE BLAT scores failed.${NORM} Aborting."
  echo "Check if file ${REDF}$TABLIST${NORM} is correct."
  echo
  exit 1
  }

# Make a new TAB file without these transcripts
echo
echo "Removing transcripts with >${CYAF}$BLATSCORE${NORM} BLAT score"
grep -v -f $TABBLAT $TAB > $TABREMOVED || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removing of transcripts with >$BLATSCORE BLAT score failed.${NORM} Aborting."
  echo "Check if files ${REDF}$TABBLAT${NORM} and ${REDF}$TAB${NORM} are correct."
  echo
  exit 1
  }
echo

{ awk '{print $1"\t"length($2)"\t"$2}' $TABREMOVED | awk '{sum+=$2}END{print sum}'; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removing of transcripts with >$BLATSCORE BLAT score failed.${NORM} Aborting."
  echo
  exit 1
  }

# Convert this TAB file to FASTA and remove the sequences with 'n'
echo
echo "Converting TAB to FASTA and removing sequences with \"n\""
{ grep -v n $TABREMOVED | sed 's/^/>/' | sed 's/\t/\n/' > $FINALA; } || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removing of transcripts with >$BLATSCORE BLAT score failed.${NORM} Aborting."
  echo "Check if file ${REDF}$TABREMOVED${NORM} is correct."
  echo
  exit 1
  }

grep -v n $TABREMOVED | awk '{print $1"\t"length($2)}' | awk '{s+=$2;a++}END{print s}' || {
  echo
  echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removing of transcripts with >$BLATSCORE BLAT score failed.${NORM} Aborting."
  echo "Check if file ${REDF}$TABREMOVED${NORM} is correct."
  echo
  exit 1
  }

# Remove unneeded temporal files - keep only *.pslx, *.fasta and *.bam
echo
echo "Removing unneeded temporal files"
if [ -n "$REFERENCEMT" ]; then
  rm $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $REFERENCEMT2* $BOWTIE2MT $FLASHOUT.extendedFrags.fastq $TAB $TABLIST $TABBLAT $TABREMOVED || {
    echo
    echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removal of temporal files failed.${NORM} Remove following files manually:"
    echo "\"$UNIQUELIST\", \"$INPUTTAB\", \"$SORTEDINPUT\","
    echo "\"$JOINEDTS\", \"$JOINEDTABS\", \"$REFERENCECP2*\","
    echo "\"$BOWTIE2CP\", \"$REFERENCEMT2*\", \"$BOWTIE2MT\","
    echo "\"$FLASHOUT.extendedFrags.fastq\", \"$TAB\", \"$TABLIST\","
    echo "\"$TABBLAT\" and \"$TABREMOVED\"."
    confirmgo
    }
  else
    rm $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $FLASHOUT.extendedFrags.fastq $TAB $TABLIST $TABBLAT $TABREMOVED || {
      echo
      echo "${REDF}${BOLD}Error!${NORM} ${CYAF}Removal of temporal files failed.${NORM} Remove following files manually:"
      echo "\"$UNIQUELIST\", \"$INPUTTAB\", \"$SORTEDINPUT\","
      echo "\"$JOINEDTS\", \"$JOINEDTABS\", \"$REFERENCECP2*\","
      echo "\"$BOWTIE2CP\", \"$FLASHOUT.extendedFrags.fastq\", \"$TAB\","
      echo "\"$TABLIST\", \"$TABBLAT\" and \"$TABREMOVED\"."
      confirmgo
      }
  fi

# List kept files which user can use for another analysis
echo
echo "${CYAF}Following files are kept for possible later usage (see manual for details):${NORM}"
echo "${BLUF}================================================================================${NORM}"
if [ -n "$REFERENCEMT" ]; then
  if [ "$RENAMEDFASTASEQENCES" == "YES" ]; then
    echo "${CYAF}0)${NORM}  List of old names of FASTA sequences in ${REDF}$INPUTFILE0${NORM} and"
    echo "    renamed FASTA sequences in ${REDF}$INPUTFILE${NORM}:"
    echo "${REDF}$TRANSCRIPTOMEFASTANAMES${NORM}"
    fi
  echo "${CYAF}1)${NORM}  Output of BLAT (removal of transcripts sharing ≥90% sequence similarity):"
  echo "${REDF}$BLATOUT${NORM}"
  echo "${CYAF}2)${NORM}  Unique transcripts in FASTA format:"
  echo "${REDF}$JOINEDFA${NORM}"
  echo "${CYAF}3)${NORM}  SAM converted to BAM (removal of reads of plastid origin):"
  echo "${REDF}$CPBAM${NORM}"
  echo "${CYAF}4)${NORM}  Genome skim data without cpDNA reads:"
  ls $FASTQNOCP*
  echo "${CYAF}5)${NORM}  SAM converted to BAM (removal of reads of mitochondrial origin):"
  echo "${REDF}$MTBAM${NORM}"
  echo "${CYAF}6)${NORM}  Genome skim data without mtDNA reads:"
  ls $FASTQNOMT*
  echo "${CYAF}7)${NORM}  Combined paired-end genome skim reads:"
  echo "${REDF}$FLASHOUT.extendedFrags.fa${NORM}"
  echo "${CYAF}8)${NORM}  Output of BLAT (matching of the unique transcripts and the filtered,"
  echo "    combined genome skim reads sharing ≥85% sequence similarity):"
  echo "${REDF}$BLATOUTFIN${NORM}"
  echo "${CYAF}9)${NORM}  Matching sequences in FASTA:"
  echo "${REDF}$BLATOUTFIN2${NORM}"
  echo "${CYAF}10)${NORM} Final FASTA sequences for usage in Geneious:"
  echo "${REDF}$FINALA${NORM}"
  else
    if [ "$RENAMEDFASTASEQENCES" == "YES" ]; then
      echo "${CYAF}0)${NORM}  List of old names of FASTA sequences in ${REDF}$INPUTFILE0${NORM} and"
      echo "    renamed FASTA sequences in ${REDF}$INPUTFILE${NORM}:"
      echo "${REDF}$TRANSCRIPTOMEFASTANAMES${NORM}"
      fi
    echo "${CYAF}1)${NORM} Output of BLAT (removal of transcripts sharing ≥90% sequence similarity):"
    echo "${REDF}$BLATOUT${NORM}"
    echo "${CYAF}2)${NORM} Unique transcripts in FASTA format:"
    echo "${REDF}$JOINEDFA${NORM}"
    echo "${CYAF}3)${NORM} SAM converted to BAM (removal of reads of plastid origin):"
    echo "${REDF}$CPBAM${NORM}"
    echo "${CYAF}4)${NORM} Genome skim data without cpDNA reads:"
    ls $FASTQNOCP*
    echo "${CYAF}5)${NORM} Combined paired-end genome skim reads:"
    echo "${REDF}$FLASHOUT.extendedFrags.fa${NORM}"
    echo "${CYAF}6)${NORM} Output of BLAT (matching of the unique transcripts and the filtered,"
    echo "   combined genome skim reads sharing ≥85% sequence similarity):"
    echo "${REDF}$BLATOUTFIN${NORM}"
    echo "${CYAF}7)${NORM} Matching sequences in FASTA:"
    echo "${REDF}$BLATOUTFIN2${NORM}"
    echo "${CYAF}8)${NORM} Final FASTA sequences for usage in Geneious:"
    echo "${REDF}$FINALA${NORM}"
  fi
echo "${BLUF}================================================================================${NORM}"
confirmgo

echo "${REDF}${BOLD}Success!${NORM}"
echo
echo "${BLUF}================================================================================${NORM}"
echo "${CYAF}Resulting FASTA was saved as${NORM}"
echo "${REDF}${BOLD}$FINALA${NORM}"
echo "for usage in Geneious (step 7 of the pipeline)."
echo "${CYAF}Use this file in next step of the pipeline. See ${REDF}README${CYAF} and ${REDF}manual${CYAF} for details.${NORM}"
echo "${BLUF}================================================================================${NORM}"
confirmgo

echo "${CYAF}Run Geneious${NORM} (tested with versions 6, 7 and 8; see ${REDF}README${NORM} and ${REDF}manual${NORM} for details):"
echo "${BLUF}================================================================================${NORM}"
echo "Import output file \"${REDF}$FINALA${NORM}\" (${CYAF}File ${REDF}|${CYAF} Import ${REDF}|${CYAF} From File...${NORM})."
echo "Select the file and go to menu ${CYAF}Tools ${REDF}|${CYAF} Align / Assemble ${REDF}|${CYAF} De Novo Assemble${NORM}."
echo "In \"${REDF}Data${NORM}\" frame select \"${CYAF}Assembe by 1st (...) Underscore${NORM}\"."
echo "In \"${REDF}Method${NORM}\" frame select Geneious Assembler (if you don't have other assemblers,"	
echo "  this option might be missing) and \"${CYAF}Medium Sensitivity / Fast${NORM}\" Sensitivity."
echo "In \"${REDF}Results${NORM}\" frame check \"${CYAF}Save assembly report${NORM}\", \"${CYAF}Save list of unused reads${NORM}\","
echo "  \"${CYAF}Save in sub-folder${NORM}\", \"${CYAF}Save contigs${NORM}\" (do not check \"${CYAF}Maximum${NORM}\") and"
echo "  \"${CYAF}Save consensus sequences${NORM}\"."
echo "${CYAF}Do not trim. Otherwise keep defaults.${NORM} ${REDF}Run it.${NORM}"
echo "Geneious may warn about possible hanging because of big file."
echo "Do not use Geneious for other tasks during the assembly."
echo
echo "${CYAF}Select all resulting contigs${NORM} (typically named \"${REDF}* Contig #${NORM}\") and ${CYAF}export them${NORM}"
echo "  (${CYAF}File ${REDF}|${CYAF} Export ${REDF}|${CYAF} Selected Documents...${NORM}) as \"${CYAF}Tab-separated table values (*.tsv)${NORM}\"."
echo "${CYAF}Save following columns${NORM} (Hold ${REDF}Ctrl${NORM} key to mark more fields): \"${CYAF}# Sequences${NORM}\","
echo "  \"${CYAF}% Pairwise Identity${NORM}\", \"${CYAF}Description${NORM}\", \"${CYAF}Mean Coverage${NORM}\", \"${CYAF}Name${NORM}\" and"
echo "  \"${CYAF}Sequence Length${NORM}\"."
echo "If this option would be inaccessible for you, export all columns."
echo "${REDF}${BOLD}Warning! Do not select and export${NORM} \"${CYAF}* Consensus Sequences${NORM}\", \"${CYAF}* Unused Reads${NORM}\""
echo "or \"${CYAF}* Assembly Report${NORM}\" - only the individual  \"${CYAF}* contig #${NORM}\" files."
echo
echo "Select items \"${REDF}Consensus Sequences${NORM}\" and \"${REDF}Unused Reads${NORM}\" and ${CYAF}export both of them"
echo "  as one FASTA file.${NORM}"
echo "Go to menu ${CYAF}File ${REDF}|${CYAF} Export ${REDF}|${CYAF} Selected Documents...${NORM} and choose ${REDF}FASTA${NORM} file type."
echo
echo "${CYAF}Use exported files from Geneious as input for part B${NORM} of the script"
echo "  (${REDF}./sondovac_part_b.sh -h${NORM}; see ${REDF}README${NORM} and/or ${REDF}PDF manual${NORM} for details)."
echo "${BLUF}================================================================================${NORM}"
echo
echo "Script exited successfully..."
echo

exit
