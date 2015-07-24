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

echo "This is part A of the pipeline."
echo
echo "This part is for filtering of raw data and their preparation for assembly in Geneious."
echo "Results of Geneious assembly are processed in part B to get the final list of low-copy nuclear probe sequences."

# Default values
# Counter if not both -i and -n options are used
CHECKMODE=0
# If not specifying explicitly otherwise (using -n), running in interactive mode
STARTINI="I"
# flash -M read length of paired-end reads
FLASHM=250
# BLAT -minIdentity between the unique transcriptomes and the genome skim data
BLATIDENT=85
# Remove transcripts with >1000 BLAT scores
BLATSCORE=1000
# Create empty variables for file names
INPUTFILE=""
REFERENCECP=""
INPUTFQ1=""
INPUTFQ2=""
REFERENCEMT=""

# Parse initial arguments
while getopts "hvulrpeinf:c:m:t:q:a:y:s:" START; do
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
      echo -e "\t\tThis is the recommended default value (the script runs interactively without explicitly using option ${BOLD}-n${NORM})."
      echo -e "\t-n\tRunning in non-interactive mode. User ${BOLD}must${NORM} provide at least five input files below:"
      echo -e "\tYou can use ${BOLD}only one${NORM} of the parameters ${BOLD}-i${NORM} or ${BOLD}-n${NORM} (not both of them)"
      echo
      echo -e "\tIf options ${BOLD}-f, -c, -m, -t${NORM} and ${BOLD}-q${NORM} are used and the script is running in interactive mode, those values will be used as defaults, but may later be overwritten."
      echo
      echo -e "\tOptions required for running in non-interactive mode:"
      echo -e "\t-f\tTranscript input file in FASTA format"
      echo -e "\t-c\tPlastom reference sequence in FASTA format"
      echo -e "\t-m\tChondrion reference sequence in FASTA format"
      echo -e "\t-t\tGenome skimming paired-end input data in FASTQ format (first file)"
      echo -e "\t-q\tGenome skimming paired-end input data in FASTQ format (second file)"
      echo -e "\tNon-interactive mode is for the situation when user the is ${BOLD}sure${NORM} that all required software is installed (see README and INSTALL for details). It is recommended for experienced users, for running on remote servers and so on."
      echo
      echo -e "\tOther optional arguments (if not provided, default values are used):"
      echo -e "\t-a\tRead length of paired-end reads (parameter -M of FLASH, see its manual for details)"
      echo -e "\t\tDefault value: 250 (allowed values are 125, 150, 250 or 300)"
      echo -e "\t-y\tBLAT score for identity between unique transcript and genome skimming data (parameter -minIdentity, check the BLAT manual for details)"
      echo -e "\t\tDefault value: 85 (range from 1 to 100; if lower target enrichment might be much less efficient)"
      echo -e "\t-s\tBLAT score (range from 100 to 100000; transcripts with higher score will be removed - very likely that these are repetitive elements)"
      echo -e "\t\tDefault value: 1000"
      echo -e "\t${BOLD}WARNING!${NORM} If parameters ${BOLD}-a, -y${NORM} or ${BOLD}-s${NORM} are not provided, default values are taken and it is not possible to change them later (not even in interactive mode)."
      echo
      echo "Examples:"
      echo "$0 -i # Basic and the most simple usage"
      echo "$0 -i -f input.fa -t reads1.fastq -q reads2.fastq"
      echo "$0 -n -f input.fa -c referencecp.fa -m referencemt.fa -t reads1.fastq -q reads2.fastq"
      echo "$0 -i -a 220 # Modify parameter ${BOLD}-a${NORM}, otherwise run interactively"
      echo "$0 -n -f input.fa -c referencecp.fa -m referencemt.fa -t reads1.fastq -q reads2.fastq -y 90"
      echo "$0 -s 950 # Note interactive mode ${BOLD}-i${NORM} is implicit and does not need to be specified explicitly"
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
      INPUTFILE=$OPTARG
      echo "Input file: $INPUTFILE"
      ;;
    c)
      REFERENCECP=$OPTARG
      echo "Plastom reference: $REFERENCECP"
      ;;
    m)
      REFERENCEMT=$OPTARG
      echo "Chondrion reference: $REFERENCEMT"
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
      case "$FLASHM" in
	125) FLASHM=125;;
	150) FLASHM=150;;
	250) FLASHM=250;;
	300) FLASHM=300;;
	*) echo
	  echo "Error! For parameter \"-a\" you did not provide any of the allowed values 125, 150, 250, or 300!"
	  echo
	  exit 1
	  ;;
      esac
      echo "Read length of paired-end reads: $FLASHM"
      ;;
    y)
      BLATIDENT=$OPTARG
      # Check if provided value makes sense
      if [[ "$BLATIDENT" =~ ^[0-9]+$ ]] && [ "$BLATIDENT" -ge 1 -a "$BLATIDENT" -le 100 ]; then
	echo "BLAT score for identity between unique transcripts and genome skimming data: $BLATIDENT"
	else
	  echo
	  echo "Error! For parameter \"-y\" you did not provide an integer of range from 1 to 100!"
	  echo
	  exit 1
	fi
      ;;
    s)
      BLATSCORE=$OPTARG
      # Check if provided value makes sense
      if [[ "$BLATSCORE" =~ ^[0-9]+$ ]] && [ "$BLATSCORE" -ge 100 -a "$BLATSCORE" -le 100000 ]; then
	echo "BLAT score: $BLATSCORE"
	else
	  echo
	  echo "Error! For parameter \"-s\" you did not provide an integer of range from 100 to 100000!"
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
  echo "${BOLD}Error!${NORM} You provided both parameters ${BOLD}-i${NORM} (interactive mode) and ${BOLD}-n${NORM} (non-interactive mode). You ${BOLD}may${NORM} use ${BOLD}only one${NORM} of them (${BOLD}either -i${NORM} or ${BOLD}-n${NORM})!"
  echo
  exit 1
  fi

# Ensure user reads introductory information
confirmgo

# Warn user this is not yet for daily usage
devrelease

# Check operating system
oscheck

# Set variables for working directory and PATH
workdirpath

# Check availability of all needed binaries

# Function to compile BLAT
function compileblat {
  {
  echo
  checktools curl &&
  checktools unzip &&
  checktools make &&
  checktools gcc &&
  echo "Downloading BLAT source code" &&
  if [ -z $MACHTYPE ]; then
    echo
    echo "Error! Variable \$MACHTYPE required by BLAT is missing. Trying to create it. If it fails, create global variable \$MACHTYPE manually, download binary or compile it on another comparable machine."
    MACHTYPE=$HOME/bin/$OSB
    else
      echo
      echo "BLAT binaries will be available in $HOME/bin/$MACHTYPE and $BIN (consider adding them into the PATH)."
    fi &&
  echo &&
  curl -o blatSrc.zip https://users.soe.ucsc.edu/~kent/src/blatSrc.zip &&
  unzip -nq blatSrc.zip &&
  cd blatSrc &&
  { mkdir -p $HOME/bin/$MACHTYPE || { echo && echo "Error! Can not create directory \"$HOME/bin/$MACHTYPE\" required by BLAT. Aborting." && echo && exit 1; }; } &&
  PATH=$PATH:$HOME/bin/$MACHTYPE &&
  { mkdir -p lib/$MACHTYPE || { echo && echo "Error! Can not create directory \"$(pwd)/lib/$MACHTYPE\" required by BLAT. Aborting." && echo && exit 1; }; } &&
  echo "Compiling BLAT from source code" &&
  make -s &&
  cd $WORKDIR &&
  cp -p $HOME/bin/$MACHTYPE/* $BIN/ &&
  echo "\"BLAT\" is available. OK"
  } || { echo "Compilation failed. Please go to https://users.soe.ucsc.edu/~kent/src/, download latest blatSrc*.zip, compile it and ensure it is in PATH" && exit 1; }
  }

# Check if BLAT is available
{ command -v blat >/dev/null 2>&1 && echo "\"BLAT\" is available. OK."; } || {
  echo
  echo "BLAT is required but not installed or available in PATH."
  echo
  if [ "$STARTINI" == "I" ]; then
    echo "Type \"S\" to compile BLAT from source available on https://users.soe.ucsc.edu/~kent/src/ (BLAT license does not allow redistributions, required if BLAT is not available for your system). Together with standard compilation tools BLAT requires libpng developmental files."
    echo "Type \"D\" to download BLAT from http://hgdownload.cse.ucsc.edu/admin/exe/ automatically for your OS (BLAT license does not allow redistributions, available for 64 bit Linux and Mac OS X, recommended)."
    echo "Type \"M\" for manual installation - script will exit and you will have to install BLAT yourselves. Check http://genome.ucsc.edu/FAQ/FAQblat.html for more information."
    read BLAT
    while :
    do
      case "$BLAT" in
	S|s)
	  compileblat
	  break
	  ;;
	D|d) 
	  if [ "$OS" == "Mac" ]; then
	    {
	    checktools curl &&
	    echo "Downloading blat binary for $OS" &&
	    curl -o blat http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat &&
	    chmod +x blat &&
	    mv blat $BIN/ &&
	    echo "\"BLAT\" is available. OK"
	    } || { echo && echo "Download of BLAT failed. Please, go to http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/blat/ and download blat binary yourselves." && echo && exit 1; }
	    break
	    elif [[ "$OS" == "Linux" && "$OSB" == "64b" ]]; then
	      {
	      echo "Downloading blat binary for $OS" &&
	      curl -o blat http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat &&
	      chmod +x blat &&
	      mv blat $BIN/ &&
	      echo "\"BLAT\" is available. OK."
	      } || { echo && echo "Download of BLAT failed. Please, go to http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/ and download blat binary yourselves." && echo && exit 1; }
	      break
	    elif [[ "$OS" == "Linux" && "$OSB" == "32b" ]]; then
	      echo "BLAT binary is not provided for 32bit Linux."
	      compileblat
	    else
	      echo "Unknown OS or OS without BLAT binary available."
	      compileblat
	    fi
	    break
	  ;;
	M|m) echo "Please, go to http://genome.ucsc.edu/FAQ/FAQblat.html and download and install BLAT and ensure it is in PATH." && echo && exit 2;;
	*) echo "Wrong option. Use S, D or M." && read BLAT;;
      esac
    done
    else
      exit 1
  fi
}

# Check if cut is available
checktools cut

# Check if awk is available
checktools awk

# Check if join is available
checktools join

# Check if sed is available
checktools sed

# Function to compile Bowtie2
function compilebowtie {
  {
  echo &&
  checktools make &&
  checktools g++ &&
  echo &&
  echo "Compiling Bowtie2 from source code" &&
  echo &&
  cd $1 &&
  make -s &&
  cp bowtie2* $BIN/ &&
  cd $WORKDIR &&
  echo &&
  echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\" and \"bowtie2-build-s\" are available. OK."
  } || { echo && echo "Compilation failed. Please go to from http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/, download bowtie2-*-source.zip, compile it and ensure it is in PATH" && echo && exit 1; }
  }

# Check if bowtie2 is available
{ command -v bowtie2 >/dev/null 2>&1 && command -v bowtie2-align-l >/dev/null 2>&1 && command -v bowtie2-align-s >/dev/null 2>&1 && command -v bowtie2-build >/dev/null 2>&1 && command -v bowtie2-build-l >/dev/null 2>&1 && command -v bowtie2-build-s >/dev/null 2>&1 && echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\" and \"bowtie2-build-s\" are available. OK."; } || {
  echo
  echo >&2 "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\" and \"bowtie2-build-s\" are required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"C\" to compile Bowtie2-2.2.5 from source available together with this script"
    echo "Type \"S\" to compile Bowtie2-2.2.5 from source code downloaded from http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/"
    echo "Type \"D\" to download Bowtie2-2.2.5 binary from http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/ for your OS"
    echo "Type \"B\" to copy Bowtie2-2.2.2.5 binary available together with the script (available for Linux and Mac OS X)"
    echo "Type \"M\" for manual installation - script will exit and you will have to install Bowtie2 yourselves"
    read BOWTIE
    while :
    do
      case "$BOWTIE" in
	C|c)
	  compilebowtie $SCRIPTDIR/src/bowtie2-2.2.5
	  break
	  ;;
	D|d)
	  if [ "$OSB" == "32b" ]; then
	    echo
	    checktools curl
	    checktools unzip
	    echo
	    echo "Bowtie binary for 32bit CPU is not available."
	    compilebowtie $SCRIPTDIR/src/bowtie2-2.2.5
	    elif [ "$OS" == "Linux" ]; then
	      {
	      echo "Downloading Bowtie2 binaries for $OS $OSB" &&
	      curl -L -o bowtie2-2.2.5-linux-x86_64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip &&
	      unzip -nq bowtie2-2.2.5-linux-x86_64.zip &&
	      cp bowtie2-2.2.5/bowtie2* $BIN/ &&
	      echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\" and \"bowtie2-build-s\" are available. OK."
	      } || { echo && echo "Download failed. Please, go to http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ and download it manually" && echo && exit 1; }
	    elif [ "$OS" == "Mac" ]; then
	      {
	      echo "Downloading Bowtie2 binaries for $OS $OSB" &&
	      curl -L -o bowtie2-2.2.5-macos-x86_64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.5/bowtie2-2.2.5-macos-x86_64.zip &&
	      unzip -nq bowtie2-2.2.5-macos-x86_64.zip &&
	      cp bowtie2-2.2.5/bowtie2* $BIN/ &&
	      echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\" and \"bowtie2-build-s\" are available. OK."
	      } || { echo && echo "Download failed. Please, go to http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ and download it manually" && echo && exit 1; }
	    elif [ "$OS" == "Windows" ]; then
	      {
	      echo "Downloading Bowtie2 binaries for $OS $OSB" &&
	      curl -L -o bowtie2-2.2.5-mingw-win64.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.5/bowtie2-2.2.5-mingw-win64.zip &&
	      unzip -nq bowtie2-2.2.5-mingw-win64.zip &&
	      cp bowtie2-2.2.5/bowtie2* $BIN/ &&
	      echo "\"bowtie2\", \"bowtie-align-l\", \"bowtie-align-s\", \"bowtie2-build\", \"bowtie2-build-l\" and \"bowtie2-build-s\" are available. OK."
	      } || { echo && echo "Download failed. Please, go to http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ and download it manually" && echo && exit 1; }
	    else
	      echo "Unknown OS or OS without binary available"
	      compilebowtie $SCRIPTDIR/src/bowtie2-2.2.5
	    fi
	  break
	  ;;
	S|s)
	  echo
	  checktools curl
	  checktools unzip
	  echo
	  echo "Downloading Bowtie2 source code"
	  curl -L -o bowtie2-2.2.5-source.zip http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.5/bowtie2-2.2.5-source.zip || { echo && echo "Download failed. Please, go to http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ and download and install it manually" && echo && exit 1; }
	  unzip -nq $WORKDIR/bowtie2-2.2.5-source.zip
	  compilebowtie bowtie2-2.2.5
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying Bowtie binaries"
	      ;;
	    Linux)
	      echo "Copying Bowtie binaries"
	      if [ "$OSB" == "64b" ]; then
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/bowtie2* $BIN/
		mkdir -p $WORKDIR/bin/share
		cp -pr $SCRIPTDIR/pkgs/linux64b/share/man $WORKDIR/bin/share/
	      else
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/bowtie2* $BIN/
		mkdir -p $WORKDIR/bin/share
		cp -pr $SCRIPTDIR/pkgs/linux32b/share/man $WORKDIR/bin/share/
	      fi
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB. Going to compile it from source code."
	      echo
	      compilebowtie $SCRIPTDIR/src/bowtie2-2.2.5
	      ;;
	  esac
	  break
	  ;;
	M|m) echo "Please, go to http://bowtie-bio.sourceforge.net/bowtie2/index.shtml and install Bowtie2 and ensure it is in PATH" && exit 2;;
	*) echo "Wrong option. Use C, D, S, B or M." && read BOWTIE;;
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
  echo "\"samtools\" is available. OK."
  }

# Check if samtools is available
{ command -v samtools >/dev/null 2>&1 && echo "\"samtools\" is available. OK."; } || {
  echo
  echo >&2 "\"samtools\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"C\" to compile Samtools-1.2 from source available together with this script. Makefile was modified not to require GNU ncurses library."
    echo "Type \"S\" to download latest developmental Samtools source from https://github.com/samtools/samtools/ and compile it. Compilation requires GNU ncurses library and is recommended only for advanced users."
    echo -e "\tIf compilation fails, check Samtools' INSTALL file for details and adjust its Makefile."
    echo "Type \"B\" to copy Samtools-1.2 binary available together with the script (available for Linux and Mac OS X)"
    echo "Type \"M\" for manual installation - script will exit and you will have to install Samtools yourselves."
    read SAMTOOLS
    while :
    do
      case "$SAMTOOLS" in
	C|c)
	  compilesamtools || { echo && echo "Compilation failed. Please go to http://www.htslib.org/download/, download samtools, compile it and ensure it is in PATH" && echo && exit 1; }
	  break
	  ;;
	S|s)
	  {
	  echo
	  checktools curl &&
	  checktools unzip &&
	  checktools make &&
	  checktools gcc &&
	  echo &&
	  echo "Downloading Samtools sources" &&
	  curl -L -o develop.zip https://github.com/samtools/samtools/archive/develop.zip &&
	  unzip -nq develop.zip &&
	  cd samtools-develop &&
	  curl -L -o develop.zip https://github.com/samtools/htslib/archive/develop.zip &&
	  unzip -nq develop.zip &&
	  echo "Compiling Samtools" &&
	  make -s HTSDIR=`pwd`/htslib-develop &&
	  make -s prefix=$WORKDIR/bin install &&
	  cd $WORKDIR &&
	  echo "\"samtools\" is available. OK."
	  } || { echo && echo "Compilation failed. Please go to http://www.htslib.org/download/, download samtools, compile it and ensure it is in PATH" && echo && exit 1; }
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying Samtools binaries"
	      ;;
	    Linux)
	      echo "Copying Samtools binaries"
	      if [ "$OSB" == "64b" ]; then
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/ace2sam $BIN/
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/blast2sam.pl $BIN/
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
	      else
	      	cp -p $SCRIPTDIR/pkgs/linux32b/bin/ace2sam $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/blast2sam.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/export2sam.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/interpolate_sam.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/maq2sam-long $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/maq2sam-short $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/md5fa $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/md5sum-lite $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/novo2sam.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/plot-bamstats $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/psl2sam.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/samtools $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/samtools.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/sam2vcf.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/seq_cache_populate.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/soap2sam.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/varfilter.py $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/wgsim $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/wgsim_eval.pl $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/zoom2sam.pl $BIN/
	      fi
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB. Going to compile it from source code."
	      echo
	      compilesamtools
	      ;;
	  esac
	  break
	  ;;
	M|m) echo "Please, go to http://www.htslib.org/ and install Samtools and ensure it is in PATH" && exit 2;;
	*) echo "Wrong option. Use C, S, B or M." && read SAMTOOLS;;
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
  echo "Compiling bam2fastq" &&
  cd $1 &&
  make -s &&
  cp bam2fastq $BIN/ &&
  cd $WORKDIR &&
  echo "\"bam2fastq\" is available. OK."
  }

# Check if bam2fastq is available
{ command -v bam2fastq >/dev/null 2>&1 && echo "\"bam2fastq\" is available. OK."; } || {
  echo
  echo >&2 "\"bam2fastq\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"C\" to compile bam2fastq from source available together with this script"
    echo "Type \"S\" to download bam2fastq source from http://gsl.hudsonalpha.org/information/software/bam2fastq and compile it"
    echo "Type \"B\" to copy bam2fastq binary available together with the script (available for Linux and Mac OS X)"
    echo "Type \"M\" for manual installation - script will exit and you will have to install bam2fastq yourselves"
    read BAM2FASTQ
    while :
    do
      case "$BAM2FASTQ" in
	C|c)
	  compilebam2fastq $SCRIPTDIR/src/bam2fastq-1.1.0 || { echo && echo "Compilation failed. Please go to http://gsl.hudsonalpha.org/information/software/bam2fastq, download bam2fastq, compile it and ensure it is in PATH" && echo && exit 1; }
	  break
	  ;;
	S|s)
	  {
	  echo &&
	  checktools curl &&
	  checktools tar &&
	  checktools gunzip &&
	  checktools make &&
	  echo &&
	  echo "Downloading bam2fastq source code" &&
	  curl -o bam2fastq-1.1.0.tgz http://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz &&
	  tar xzvf bam2fastq-1.1.0.tgz &&
	  compilebam2fastq bam2fastq-1.1.0/
	  } || { echo && echo "Compilation failed. Please go to http://gsl.hudsonalpha.org/information/software/bam2fastq, download bam2fastq, compile it and ensure it is in PATH" && echo && exit 1; }
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying bam2fastq binary"
	      ;;
	    Linux)
	      echo "Copying bam2fastq binary"
	      if [ "$OSB" == "64b" ]; then
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/bam2fastq $BIN/
	      else
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/bam2fastq $BIN/
	      fi
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB. Going to compile it from source code."
	      echo
	      compilebam2fastq $SCRIPTDIR/src/bam2fastq-1.1.0 || { echo && echo "Compilation failed. Please go to http://gsl.hudsonalpha.org/information/software/bam2fastq, download bam2fastq, compile it and ensure it is in PATH" && echo && exit 1; }
	      ;;
	  esac
	  break
	  ;;
	M|m) echo "Please, go to http://www.htslib.org/ and install Samtools and ensure it is in PATH" && exit 2;;
	*) echo "Wrong option. Use C, S, B or M." && read BAM2FASTQ;;
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
  echo "Compiling FLASH from source code" &&
  echo &&
  cd $1 &&
  make -s &&
  cp flash $BIN/ &&
  cd $WORKDIR &&
  echo &&
  echo "\"FLASH\" is available. OK."
  } || { echo && echo "Compilation failed. Please go to http://sourceforge.net/projects/flashpage/files/, download FLASH-*.tar.gz, compile it and ensure it is in PATH" && echo && exit 1; }
  }

# Check if FLASH is available
{ command -v flash >/dev/null 2>&1 && echo "\"FLASH\" is available. OK."; } || {
  echo
  echo >&2 "\"FLASH\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"C\" to compile FLASH from source available together with this script"
    echo "Type \"S\" to download FLASH source from http://sourceforge.net/projects/flashpage/ and compile it"
    echo "Type \"D\" to download FLASH binary from http://sourceforge.net/projects/flashpage/ (available only for Windows)"
    echo "Type \"B\" to copy FLASH 1.2.11 binary available together with the script (available for Linux and Mac OS X)"
    echo "Type \"M\" for manual installation - script will exit and you will have to install FLASH yourselves"
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
	  checktools curl
	  checktools tar
	  checktools gunzip
	  echo
	  echo "Downloading FLASH source code"
	  curl -L -o FLASH-1.2.11.tar.gz http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz || { echo && echo "Download failed. Please go to http://sourceforge.net/projects/flashpage/files/, download FLASH-*.tar.gz, compile it and ensure it is in PATH" && echo && exit 1; }
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
	      checktools curl
	      checktools unzip
	      echo
	      echo "Downloading FLASH for $OS"
	      curl -L -o FLASH-1.2.11-windows-bin.zip http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11-windows-bin.zip || { echo && echo "Download failed. Please go to http://sourceforge.net/projects/flashpage/files/, download FLASH-*windows-bin.zip, unpack it and ensure it is in PATH" && echo && exit 1; }
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
	      ;;
	    Linux)
	      echo "Copying FLASH binary"
	      if [ "$OSB" == "64b" ]; then
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/flash $BIN/
	      else
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/flash $BIN/
	      fi
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB. Going to compile it from source code."
	      echo
	      compileflash $SCRIPTDIR/src/FLASH-1.2.11
	      ;;
	  esac
	  break
	  ;;
	M|m) echo "Please, go to http://ccb.jhu.edu/software/FLASH/ and install FLASH and ensure it is in PATH" && exit 2;;
	*) echo "Wrong option. Use C, S, D, B or M." && read FLASH;;
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
  checktools pkg-config &&
  echo &&
  echo "Compiling FASTX-Toolkit from source code" &&
  echo &&
  cd $1 &&
  ./configure --prefix=$WORKDIR/bin &&
  make -s &&
  make -s check &&
  make -s install &&
  make -s installcheck &&
  if [ "$OSB" == "64b" ]; then
    PKGTEXTUTILIB="lib64"
    else
      PKGTEXTUTILIB="lib"
      fi
  export PKG_CONFIG_PATH=$WORKDIR/bin/$PKGTEXTUTILIB/pkgconfig &&
  cd $2 &&
  ./configure --prefix=$WORKDIR/bin &&
  make -s &&
  make -s check &&
  make -s install &&
  make -s installcheck &&
  cd $WORKDIR &&
  echo &&
  echo "\"fastq_to_fasta\" is available. OK."
  } || { echo && echo "Compilation failed. Please go to http://hannonlab.cshl.edu/fastx_toolkit/download.html, download libgtextutils-*.tar.gz and fastx_toolkit-*.tar.bz2, compile them and ensure they are in PATH" && echo && exit 1; }
  }

# Check if fastq_to_fasta is available
{ command -v fastq_to_fasta >/dev/null 2>&1 && echo "\"fastq_to_fasta\" is available. OK."; } || {
  echo
  echo >&2 "\"fastq_to_fasta\" is required but not installed or available in PATH."
  if [ "$STARTINI" == "I" ]; then
    echo
    echo "Type \"C\" to compile FASTX-Toolkit from source available together with this script"
    echo "Type \"S\" to download FASTX-Toolkit source from http://hannonlab.cshl.edu/fastx_toolkit/ and compile it"
    echo "Type \"B\" to FASTX-Toolkit 0.0.14 binary available together with the script (available for Linux and Mac OS X)"
    echo "Type \"M\" for manual installation - script will exit and you will have to install FASTX-Toolkit yourselves"
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
	  checktools curl
	  checktools tar
	  checktools gunzip
	  checktools bunzip2
	  echo
	  echo "Downloading FASTX source code"
	  echo
	  { curl -L -o libgtextutils-0.7.tar.gz https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz &&
	  curl -L -o fastx_toolkit-0.0.14.tar.bz2 https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2; } || { echo && echo "Download failed. Please, go to https://github.com/agordon/fastx_toolkit/releases/download/ and download libgtextutils and fastx_toolkit, install them and ensure they are in PATH" && echo && exit 1; }
	  tar xzf libgtextutils-0.7.tar.gz
	  tar xjvf fastx_toolkit-0.0.14.tar.bz2
	  compilefastx $WORKDIR/libgtextutils-0.7/ $WORKDIR/fastx_toolkit-0.0.14/
	  break
	  ;;
	B|b)
	  case "$OS" in
	    Mac)
	      echo "Copying FASTX binaries"
	      ;;
	    Linux)
	      echo "Copying FASTX binaries"
	      if [ "$OSB" == "64b" ]; then
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/fasta_* $BIN/
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/fastq_* $BIN/
		cp -p $SCRIPTDIR/pkgs/linux64b/bin/fastx_* $BIN/
		mkdir -p $WORKDIR/bin/include $WORKDIR/bin/lib64 $WORKDIR/bin/share
		cp -pr $SCRIPTDIR/pkgs/linux64b/include/gtextutils $WORKDIR/bin/include/
		cp -pr $SCRIPTDIR/pkgs/linux64b/lib64/pkgconfig $WORKDIR/bin/lib64/
		cp -p $SCRIPTDIR/pkgs/linux64b/lib64/libgtextutils* $WORKDIR/bin/lib64/
		cp -pr $SCRIPTDIR/pkgs/linux64b/share/aclocal $WORKDIR/bin/share/
	      else
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/fasta_* $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/fastq_* $BIN/
		cp -p $SCRIPTDIR/pkgs/linux32b/bin/fastx_* $BIN/
		mkdir -p $WORKDIR/bin/include $WORKDIR/bin/lib $WORKDIR/bin/share
		cp -pr $SCRIPTDIR/pkgs/linux32b/include/gtextutils $WORKDIR/bin/include/
		cp -pr $SCRIPTDIR/pkgs/linux32b/lib/pkgconfig $WORKDIR/bin/lib/
		cp -p $SCRIPTDIR/pkgs/linux32b/lib32/libgtextutils* $WORKDIR/bin/lib32/
		cp -pr $SCRIPTDIR/pkgs/linux32b/share/aclocal $WORKDIR/bin/share/
	      fi
	      ;;
	    *) echo
	      echo "Binary is not available for $OS $OSB. Going to compile it from source code."
	      echo
	      compilefastx $SCRIPTDIR/src/libgtextutils-0.7/ $SCRIPTDIR/src/fastx_toolkit-0.0.14/
	      ;;
	  esac
	  break
	  ;;
	M|m) echo "Please go to http://hannonlab.cshl.edu/fastx_toolkit/download.html, download libgtextutils-*.tar.gz and fastx_toolkit-*.tar.bz2, compile them and ensure they are in PATH" && exit 2;;
	*) echo "Wrong option. Use C, S, B or M." && read FASTX;;
      esac
    done
  else
    exit 1
  fi
  }

# Check if grep is available
checktools grep

# Check if cat is available
checktools cat

# Check if perl is available
checktools perl

CHECKFILEREADOUT=""
readinputfile -f "FASTA transcriptomic data" $INPUTFILE
INPUTFILE=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

readinputfile -c "FASTA reference plastom data" $REFERENCECP
REFERENCECP=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

readinputfile -m "FASTQ reads data file 1" $INPUTFQ1
INPUTFQ1=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

readinputfile -t "FASTQ reads data file 2" $INPUTFQ2
INPUTFQ2=$CHECKFILEREADOUT
CHECKFILEREADOUT=""

readinputfile -q "FASTA reference mitochondrion data" $REFERENCEMT
REFERENCEMT=$CHECKFILEREADOUT
CHECKFILEREADOUT=""
echo

# Input file in FASTA format
echo "Input file: $INPUTFILE"
# Output of BLAT
BLATOUT="${INPUTFILE%.*}_blat_unique_transcripts.psl"
# List of unique transcripts - temporary file - will be deleted
UNIQUELIST="${INPUTFILE%.*}_trans-trans_unique_transcripts_sorted.txt"
# Input file converted into TXT - temporary file - will be deleted
INPUTTAB="${INPUTFILE%.*}.txt"
# Sorted input file in TXT - temporary file - will be deleted
SORTEDINPUT="${INPUTFILE%.*}_sorted.txt"
# Joined unique transcriptomes - temporary file - will be deleted
JOINEDTS="${INPUTFILE%.*}_unique_transcriptomes_trans-trans_plus_sequence.txt"
# Joined unique transcriptomes in tabular format - temporary file - will be deleted
JOINEDTABS="${JOINEDTS%.*}_tabs.txt"
# Joined unique transcriptomes in FASTA format
JOINEDFA="${INPUTFILE%.*}_unique_transcripts.fasta"
# Input - reference genome - cpDNA
echo "Input file: $REFERENCECP"
# Reference genome - plastome index - temporary file - will be deleted
REFERENCECP2="${REFERENCECP%.*}.cp"
# Input cpDNA reads in FASTQ
echo "Input file: $INPUTFQ1"
echo "Input file: $INPUTFQ2"
# cpDNA reads mapped to reference - temporary file - will be deleted
BOWTIE2CP="${INPUTFILE%.*}_genome_skim_data_no_cp_reads.sam"
# SAM converted into BAM
CPBAM="${INPUTFILE%.*}_genome_skim_data_no_cp_reads.bam"
# Removed cpDNA reads
FASTQNOCP="${INPUTFILE%.*}_genome_skim_data_no_cp_reads"
# Input - reference genome - mtDNA
echo "Input file: $REFERENCEMT"
# Reference genome - mitochondrion index - temporary file - will be deleted
REFERENCEMT2="${REFERENCEMT%.*}.mt"
# mtDNA reads mapped to reference - temporary file - will be deleted
BOWTIE2MT="${INPUTFILE%.*}_genome_skim_data_no_cp_no_mt_reads.sam"
# SAM converted into BAM
MTBAM="${INPUTFILE%.*}_genome_skim_data_no_cp_no_mt_reads.bam"
# Removed mtDNA reads
FASTQNOMT="${INPUTFILE%.*}_genome_skim_data_no_cp_no_mt_reads"
# Combined paired-end reads
FLASHOUT="${INPUTFILE%.*}_combined_reads_co_cp_no_mt_reads"
# BLAT between the unique transcriptomes and the genome skimming data
BLATOUTFIN="${INPUTFILE%.*}_blat_unique_transcripts_versus_genome_skim_data.pslx"
# Assembled sequences in contigs
BLATOUTFIN2="${INPUTFILE%.*}_blat_unique_transcripts_versus_genome_skim_data.fasta"
# FASTA converted into TABe - temporary file - will be deleted
TAB="${INPUTFILE%.*}_final.tab"
# Number of times each transcript hit a genome skim read - will be deleted
TABLIST="${INPUTFILE%.*}_transcript_hits.txt"
# Listed transcripts with >1000 BLAT score - will be deleted
TABBLAT="${INPUTFILE%.*}_1k_transcripts"
# Transcribds without >1000 BLAT score - will be deleted
TABREMOVED="${INPUTFILE%.*}_1k_transcripts-removed.tab"
# Final FASTA for usage in Geneious
FINALA="${INPUTFILE%.*}_blat_unique_transcripts_versus_genome_skim_data-no_missing_fin.fsa"

# Part 1: Obtain unique transcripts.

# BLAT between the transcriptome itself
echo
echo "Launching BLAT between the transcriptome itself"
blat -t=dna -q=dna -noHead -out=psl $INPUTFILE $INPUTFILE $BLATOUT || { echo && echo "${BOLD}Error!${NORM} BLAT failed. Aborting. Check if $INPUTFILE is correct." && echo && exit 1; }
echo

# Filtering for unique transcripts
echo
echo "Filtering for unique transcripts:"
{ cut -f 10 $BLATOUT | uniq -c | awk '{print$1}' | sort | uniq -c; } || { echo && echo "${BOLD}Error!${NORM} Filtering of BLAT output failed. Aborting. Check if files $INPUTFILE and $BLATOUT are correct." && echo && exit 1; }
echo
echo "Filtered transcripts saved for possible later usage as $BLATOUT for possible later usage"

# Make a list of these unique transcripts (names and sequences) and convert this file to .fasta
echo
echo "Making list of unique transcripts"
cut -f 10 $BLATOUT | uniq -c | awk '{if($1==1){print $0}}' | awk '{print $2}' | awk '{printf "%05d\n", $0;}' > $UNIQUELIST || { echo && echo "${BOLD}Error!${NORM} Making list of unique transcripts failed. Aborting. Check if files $BLATOUT and $INPUTFILE are correct." && echo && exit 1; }
echo

# In order to use the join command the original transcriptome file has to be converted to .txt with JOINEDFA2TAB, the transcript numbers have to be adjusted and the file sorted:
echo
echo "Converting original data into TXT for subsequent joining"
### REWRITE !!!
perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count JOINEDFA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' $INPUTFILE > $INPUTTAB || { echo && echo "${BOLD}Error!${NORM} Conversion of $INPUTFILE failed. Check if it is valid FASTA file. Aborting." && echo && exit 1; }
echo
echo "Sorting unique transcripts"
{ awk '{$1=sprintf("%05d", $1); print $0}' $INPUTTAB | sort > $SORTEDINPUT; } || { echo && echo "${BOLD}Error!${NORM} Sorting of unique transcripts failed. Aborting. Check if $INPUTFILE is correct FASTA file and check file $INPUTTAB." && echo && exit 1; }
echo
# Apply the join command
join -j1 $UNIQUELIST $SORTEDINPUT > $JOINEDTS || { echo && echo "${BOLD}Error!${NORM} Joining failed. Aborting. Check files $UNIQUELIST and $SORTEDINPUT if they have same number of lines." && echo && exit 1; }

# Convert to .fasta
echo
echo "Converting to FASTA"
sed -r 's/ /\t/g' $JOINEDTS > $JOINEDTABS || { echo && echo "${BOLD}Error!${NORM} Conversion of $JOINEDTS to FASTA failed. Aborting. Check file $JOINEDTS." && echo && exit 1; }
echo
awk '{print ">"$1"\n"$2}' $JOINEDTABS > $JOINEDFA || { echo && echo "${BOLD}Error!${NORM} Conversion of $JOINEDTS to FASTA failed. Aborting. Check file $JOINEDTABS." && echo && exit 1; }
echo "Joined transcripts written in FASTA format as $JOINEDFA for possible later usage"


# Part 2: Find genome skimming data (only nuclear reads) which align to the unique transcripts

# Get rid of the chloroplast and mitochondrial reads in the genome skimming data

# Chloroplast reads

# Create a reference plastome index with BOWTIE2-BUILD
echo
echo "Creating a reference plastome index"
echo
bowtie2-build $REFERENCECP $REFERENCECP2 || { echo && echo "${BOLD}Error!${NORM} Creating a reference plastome index with bowtie2-build failed. Aborting. Check if file $REFERENCECP is correct." && echo && exit 1; }
echo

# Map the cpDNA reads to the reference plastome with BOWTIE2
echo
echo "Mapping cpDNA reads to the reference plastome. This may take longer time."
echo
bowtie2 -x $REFERENCECP2 -1 $INPUTFQ1 -2 $INPUTFQ2 -S $BOWTIE2CP || { echo && echo "${BOLD}Error!${NORM} Mapping cpDNA reads to the reference plastome with bowtie2 failed. Aborting. Check if files $REFERENCECP, $REFERENCECP2, $INPUTFQ1 and $INPUTFQ2 are correct." && echo && exit 1; }
echo
echo "Mapping finished"

# Convert .sam to .bam with SAMTOOLS
echo
echo "Converting SAM to BAM. This may take longer time."
samtools view -bT $REFERENCECP $BOWTIE2CP > $CPBAM || { echo && echo "${BOLD}Error!${NORM} Conversion of SAM to BAM with samtools failed. Aborting. Check files $REFERENCECP and $BOWTIE2CP." && echo && exit 1; }
echo "Converted file saved as $CPBAM for possible later usage"

# Remove the cpDNA reads with BAM2FASTQ
echo
echo "Removing cpDNA reads. This may take longer time."
bam2fastq --no-aligned $CPBAM -o $FASTQNOCP#.fq || { echo && echo "${BOLD}Error!${NORM} Removal of cpDNA reads with bam2fastq failed. Aborting. Check if file $CPBAM is correct." && echo && exit 1; }
echo
echo "Removed reads saved for possible later usage as"
ls -1 $FASTQNOCP*

# Mitochondrial reads

# Create a reference mitochondrion index with BOWTIE2-BUILD
echo
echo "Creating a reference mitochondrion index"
echo
bowtie2-build $REFERENCEMT $REFERENCEMT2 || { echo && echo "${BOLD}Error!${NORM} Creating of reference mitochondrion index with bowtie2-build failed. Aborting. Check files $REFERENCEMT and $REFERENCEMT2." && echo && exit 1; }
echo

# Map the mtDNA reads to the reference mitochondrion with BOWTIE2
echo
echo "Mapping mtDNA reads to reference mitochondrion. This may take longer time."
echo
bowtie2 -x $REFERENCEMT2 -1 `echo $FASTQNOCP`_1.fq -2 `echo $FASTQNOCP`_2.fq -S $BOWTIE2MT || { echo && echo "${BOLD}Error!${NORM} Mapping mtDNA reads to reference mitochondrion with bowtie2 failed. Aborting. Check files $REFERENCEMT2, `echo $FASTQNOCP`_1.fq and `echo $FASTQNOCP`_2.fq." && echo && exit 1; }
echo

# Convert .sam to .bam with SAMTOOLS
echo
echo "Converting SAM to BAM. This may take longer time."
samtools view -bT $REFERENCEMT $BOWTIE2MT > $MTBAM || { echo && echo "${BOLD}Error!${NORM} Conversion of SAM to BAM failed. Aborting. Check files $REFERENCEMT and $BOWTIE2MT." && echo && exit 1; }
echo "Converted file saved as $MTBAM for possible later usage"

# Remove the mtDNA reads with BAM2FASTQ
echo
echo "Removing mtDNA reads. This may take longer time."
bam2fastq --no-aligned $MTBAM -o $FASTQNOMT#.fq || { echo && echo "${BOLD}Error!${NORM} Removal of mtDNA reads failed. Aborting. Check if file $MTBAM is correct." && echo && exit 1; }
echo

# Combine the paired-end reads with FLASH
echo
echo "Combining paired-end reads"
echo
flash -M $FLASHM `echo $FASTQNOMT`_1.fq `echo $FASTQNOMT`_2.fq -o $FLASHOUT || { echo && echo "${BOLD}Error!${NORM} Combining paired-end reads failed. Aborting. Check files $FLASHM, `echo $FASTQNOMT`_1.fq and `echo $FASTQNOMT`_2.fq." && echo && exit 1; }
echo

# Convert FASTQ file to FASTA
echo
echo "Converting FASTQ to FASTA. This may take longer time."
fastq_to_fasta -i $FLASHOUT.extendedFrags.fastq -o $FLASHOUT.extendedFrags.fa || { echo && echo "${BOLD}Error!${NORM} Conversion of FASTQ to FASTA failed. Aborting. Check if file $FLASHOUT.extendedFrags.fastq is correct." && echo && exit 1; }
echo "Converted FASTA saved as $FLASHOUT.extendedFrags.fa for possible later usage"

# BLAT between the unique transcripts and the genome skimming data
echo
echo "BLAT between the unique transcripts and the genome skimming data. This may take longer time."
blat -t=dna -q=dna -minIdentity=$BLATIDENT -out=pslx $JOINEDFA $FLASHOUT.extendedFrags.fa $BLATOUTFIN || { echo && echo "${BOLD}Error!${NORM} BLAT between the unique transcripts and the genome skimming data failed. Aborting. Check if files $JOINEDFA and $FLASHOUT.extendedFrags.fa are correct." && echo && exit 1; }
echo "Output saved as $BLATOUTFIN for possible later usage"


# Part 3: Assemble the obtained sequences in contigs (part A)

# This will be done in GENEIOUS. Modification of the .pslx file is needed (remove headers, select the field with the transcript (target) sequence names and the field with the query sequences, convert to .fasta)
echo
echo "Modifying PSLX BLAT output for usage in Geneious"
{ sed 1,5d $BLATOUTFIN | cut -f14,22 | awk '{n=split($2,a,",");for(i=1;i<=n;i++)print $1"_"NR"_"i,a[i]}' | sed s/^/'>'/ | sed s/' '/\\n/ > $BLATOUTFIN2; } || { echo && echo "${BOLD}Error!${NORM} Modifying PSLX BLAT output failed. Aborting. Check if file $BLATOUTFIN is correct." && echo && exit 1; }
echo
echo "Modified file saved as $BLATOUTFIN2 for possible later usage"

# Convert .fasta to .tab
echo
echo "Converting FASTA to TAB"
# REWRITE!!!
perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' $BLATOUTFIN2 > $TAB || { echo && echo "${BOLD}Error!${NORM} Conversion of FASTA to TAB failed. Aborting. Check if file $BLATOUTFIN2 is correct." && echo && exit 1; }
echo
{ awk '{print $1"\t"length($2)"\t"$2}' $TAB | awk '{sum+=$2}END{print sum}'; } || { echo && echo "${BOLD}Error!${NORM} Conversion of FASTA to TABConversion of FASTA to TAB failed. Aborting. Check if file $TAB is correct." && echo && exit 1; }

# Remove transcripts with >1000 BLAT scores (very likely that these are repetitive elements)

# Count the number of times each transcript hit a genome skim read
echo
echo "Counting number of times each transcript hit a genom skim read"
{ cut -f1 -d_ $TAB | sort | uniq -c | sort -n -r > $TABLIST; } || { echo && echo "${BOLD}Error!${NORM} Counting of number of times each transcript hit a genom skim read failed. Aborting. Check file $TAB." && echo && exit 1; }

# List of the transcripts with >1000 BLAT scores
echo
echo "Listing transcripts with >$BLATSCORE BLAT scores"
{ awk '$1>"'"$BLATSCORE"'"' $TABLIST | awk '{print $2}' > $TABBLAT; } || { echo && echo "${BOLD}Error!${NORM} Listing of transcripts with >$BLATSCORE BLAT scores failed. Aborting. Check if file $TABLIST is correct." && echo && exit 1; }

# Make a new .tab file without these transcripts
echo
echo "Removing transcripts with >$BLATSCORE BLAT score"
grep -v -f $TABBLAT $TAB > $TABREMOVED || { echo && echo "${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting. Check files $TABBLAT and $TAB" && echo && exit 1; }
echo
{ awk '{print $1"\t"length($2)"\t"$2}' $TABREMOVED | awk '{sum+=$2}END{print sum}'; } || { echo && echo "${BOLD}Error!${NORM} XXX failed. Aborting." && echo && exit 1; }

# Convert this .tab file to .fasta and remove the sequences with 'n'
echo
echo "Converting TAB to FASTA and removing sequences with \"n\""
{ grep -v n $TABREMOVED | sed s/^/'>'/ | sed s/\\t/\\n/ > $FINALA; } || { echo && echo "${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting. Check file $TABREMOVED." && echo && exit 1; }
grep -v n $TABREMOVED | awk '{print $1"\t"length($2)}' | awk '{s+=$2;a++}END{print s}' || { echo && echo "${BOLD}Error!${NORM} Removing of transcripts with >$BLATSCORE BLAT score failed. Aborting. Check file $TABREMOVED." && echo && exit 1; }

# Remove unneeded temporal files - keep only *.pslx, *.fasta, *.bam
echo "Removing unneeded temporal files"
rm $UNIQUELIST $INPUTTAB $SORTEDINPUT $JOINEDTS $JOINEDTABS $REFERENCECP2* $BOWTIE2CP $REFERENCEMT2* $BOWTIE2MT $FLASHOUT* $TAB $TABLIST $TABBLAT $TABREMOVED

echo
echo "Success!"
echo
echo "Resulting FASTA was saved as ${BOLD}$FINALA${NORM} for usage in Geneious"
echo
echo "Run Geneious (tested with versions 6, 7 and 8). See http://www.geneious.com/ for download and purchase."
echo "Import output file \"$FINALA\" (File | Import | From File...)."
echo "Select the file and go to menu Tools | Align / Assemble | De Novo Assemble."
echo "In \"Data\" frame select \"Assembe by 1st (...) Underscore\"."
echo "In \"Method\" frame select Geneious Assembler (if you don't have other assemblers, this option might be missing) and \"Medium Sensitivity / Fast\" Sensitivity."
echo "In \"Results\" frame check \"Save assembly report\", \"Save list of unused reads\", \"Save in sub-folder\", \"Save contigs\" (do not check \"Maximum\") and \"Save consensus sequences\"."
echo "Do not trim. Otherwise keep defaults. Run it."
echo "Geneious may warn about possible hanging because of big file. Do not use Geneious for other tasks during the assembly."
echo
echo "Select all resulting contigs (typically named \"* Contig #\") and export them (File | Export | Selected Documents...) as \"Tab-separated table values (*.tsv)\"."
echo "Save following columns (Hold Ctrl key to mark more fields): \"# Sequences\", \"% Pairwise Identity\", \"Description\", \"Mean Coverage\", \"Name\" and \"Sequence Length\"."
echo "If this option would be inaccessible for you, export all columns."
echo
echo "Select items \"Consensus Sequences\" and \"Unused Reads\" and export them as one FASTA."
echo "Go to menu File | Export | Selected Documents... and choose FASTA file type."
echo
echo "Use exported files from Geneious as input for part B of the script (see README for details)."
echo
echo "Script exited successfully..."
echo

exit
