Sondovač changelog
================================================================================

Sondovač is a script to create orthologous low-copy nuclear probes from transcriptome and genome skim data for target enrichment.

# Version 1.X regular release released 2021-MM-DD (WIP)

* Bug fixes, more robust coding.
* Updated documentation.
* Removed colorization of the output, simplified output for user, improved coding standard.
* Removed build-in software installation.
* Removed interactive mode, simplified code.

# Version 1.3 regular release released 2017-12-18

* bam2fastq is dropped in favour of samtools fastq. No plans to use Picard anymore (part A).
* Simplified INSTALL and README not to only copy PDF manual.
* Corrected output of part A - ensure to have always valid FASTA.
* Automatically remove putative plastid sequences from final probe set (part B), list of all probes (with putative plastid sequences) and list of putative plastid sequences are available.
* Updated software distributed with Sondovač, updated respective sections of PDF manual.
* Removed FASTX Toolkit, conversion from FASTQ to FASTA is done by simple shell function.
* Improved handling of input/output files when stored in several different directories.
* Tested with Geneious 10, improved description of Geneious usage in the PDF manual.
* Improved PDF manual.

# Version 1.2 regular release released 2016-06-28

* Fixed wrong probe design summary statistics (at the end of part B).
* Added support for Geneious 9.
* More detailed error reporting (especially in part B).
* Various smaller fixes.
* Updated accompanying software to recent versions.

# Version 1.1 regular release released 2016-03-15

* Checking if input FASTA files are interleaved or not (required) and, if needed, FASTA files are converted not to be interleaved.
* Added requirement of Python.
* Removed BAM files (part A) - they are not kept anymore.
* Language checking and enhancements of documentation.
* Modified CD-HIT-EST.
* Improved summary of the probe design in part B.
* More possibilities for minimum total locus length.
* Various smaller fixes.

# Version 1.0 regular release released 2016-01-12

* Renaming of input FASTA sequences names is required - it ensures correct working of part B.
* Added check if input files were created on Windows - if so, they are converted into UNIX style EOL.
* Various smaller fixes.
* Better showing of the information in part B.
* Enhanced documentation.

# Version 0.99 release candidate released 2015-12-08

* Fixed error with some input files for part B.
* Finished colorization of command-line user interface.
* Added possibility to set minimal exon length of the loci.
* Various fixes and UI enhancements.
* Improved documentation.

# Version 0.95 beta released 2015-11-27

* Offer the possibility to choose between transcripts or genome skim sequences for further processing.
* Colorization of command-line user interface (incomplete).
* Added possibility to change -minIdentity parameter of BLAT in step 11, part B.
* Fixed problems with some transcriptome input files.
* Added possibility to set custom bait length.
* Added information about article in MER introducing Sondovač.

# Version 0.9 beta released 2015-10-23

* Highly enhanced part B.
* Better handling of variable output from Geneious.
* Possibility to specify custom output files name.
* Full support for Linux distributions using DEB - Debian, Ubuntu, Linux Mint and derivatives.
* Enhanced documentation.
* Support for Mac OS X, package management using Homebrew.
* Support for RedHat based Linux distributions - Fedora, Centos and Scientific Linux and derivatives.
* Better compilation and installation of required software.
* For downloading automatically select whether to use wget (preferred) or curl.
* Various fixes.

# Version 0.8 alpha released 2015-10-09

* Usage of mitochondrial reference sequence is optional.
* Better formatting of script messages.
* Various fixes and enhancements.

# Version 0.7 alpha released 2015-10-06

* Fixed reported problems with sed differences among Linux and Mac OS X.
* Added more exhaustive documentation.
* Various fixes and enhancements.

# Version 0.6 alpha released 2015-08-10

* Fixed problems with some versions of output of Geneious.
* Better compilation and installation of required additional software packages.
* Various fixes and enhancements.

# Version 0.5 alpha released 2015-07-24

* First public release, early alpha stage.

