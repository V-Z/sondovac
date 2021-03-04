Sondovač install instructions
================================================================================

Sondovač is a script to create orthologous low-copy nuclear probes from transcriptome and genome skim data for target enrichment.

# Requirements to run Sondovač

Sondovač is tested on openSUSE and Debian Linux distributions. It will run on any Linux/UNIX operating system.

In order to run Sondovač you need a UNIX-based operating system (preferably Linux) equipped with BASH or compatible shell interpreter (this should by default be available for any Linux distribution. You should use the current operating system version supported by upstream. Older operating systems can have different versions of shell and system libraries, which can cause various problems and incompatibilities.

Sondovač uses several scientific software packages (namely [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [CD-HIT](http://weizhongli-lab.org/cd-hit/), [FLASH](https://sourceforge.net/projects/flashpage/), [Geneious](https://www.geneious.com/), and [SAMtools](https://www.htslib.org/)), and basic UNIX tools. Sondovač will check if those programs are installed --- available in the PATH (i.e. if the shell application can locate and launch respective binaries). If you have those packages installed (in current versions), ensure their binaries are in PATH. This should not be a problem for basic tools available in any UNIX-based operating system, as basic installation usually contains all needed tools. If you lack some of the required tools, the script will notify you, and you will have to install them manually. If this is needed, check the documentation for your operating system. Users of macOS can install those applications also using [Homebrew](https://brew.sh/). For compilation of needed software you need GNU G++, GNU GCC, libpng developmental files, and zlib developmental files. Ensure that you have those tools available --- they should be readily available for any UNIX-based operating system.

**`sondovac_part_a.sh` requires the following software packages:**

* BLAT: <https://genome.ucsc.edu/FAQ/FAQblat.html#blat3> --- easiest is to download binary from <https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/>
* Bowtie2: <https://sourceforge.net/projects/bowtie-bio/files/bowtie2/>
* SAMtools: <https://github.com/samtools/samtools/releases>
* FLASH: <https://sourceforge.net/projects/flashpage/files/>

**`sondovac_part_b.sh` requires the following software packages:**

* BLAT: <https://genome.ucsc.edu/FAQ/FAQblat.html#blat3> --- easiest is to download binary from <https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/>
* CD-HIT: <https://github.com/weizhongli/cdhit/releases>

**Geneious** is required for step 7 of the pipeline. See below, `README.md` and PDF manual for details.

Basic UNIX tools are required to run Sondovač. They are usually readily available in UNIX systems (coreutils), so there is usually no need to install them manually. macOS is missing some tools and for others (typically sed, grep or awk) contains outdated BSD versions. We recommend to install [Homebrew](https://brew.sh/) and use it to install all needed tools (coreutils, gnu-sed, gawk, grep, bash, and also the scientific software). See the PDF manual for details about tools required by Sondovač and their manual installation.

# Install Sondovač

Download the latest version of Sondovač from <https://github.com/V-Z/sondovac/releases/>, unpack it and start by

    ./sondovac_part_a.sh -h

to see basic usage instructions. See `README.md` and PDF manual for more information.

# Help for usage of terminal

If you are not familiar with the usage of command line, familiarise yourself with some basic tutorials first. You can try some of these:

* <https://soubory.trapa.cz/linuxcourse/>
* <https://doc.opensuse.org/documentation/leap/startup/html/book-opensuse-startup/part-bash.html>
* <https://help.ubuntu.com/community/UsingTheTerminal>
* <https://www.gnu.org/software/bash/manual/> (advanced --- full reference manual)
* <https://www.debian.org/doc/manuals/debian-reference/ch01.en.html>
* <https://en.wikibooks.org/wiki/Guide_to_Unix>
* <https://tldp.org/LDP/Bash-Beginners-Guide/html/Bash-Beginners-Guide.html>
* <https://linuxcourse.rutgers.edu/documents/Bash-Beginners-Guide/>
* <https://ryanstutorials.net/linuxtutorial/>
* <https://www.hypexr.org/bash_tutorial.php>
* <https://mywiki.wooledge.org/BashGuide>

# Steps of the pipeline

Sondovač workflow is divided into three parts (see `README.md` and PDF manual for details):

- **1)** Raw input data are analyzed by `sondovac_part_a.sh`.
- **2)** Sequences obtained in part A are assembled by [Geneious](https://www.geneious.com/) in a separate step by the user.
- **3)** Final probes are produced by `sondovac_part_b.sh`.

For part (2) of the script the user must have [Geneious](https://www.geneious.com/). We plan to replace it with a free open-source command line tool in a future release of Sondovač. Visit <https://www.geneious.com/> for download, purchase, installation and usage of Geneious.

