Sondovač Basic help
================================================================================

**Sondovač is a script to create orthologous low-copy nuclear probes from transcriptome and genome skim data for target enrichment.**

# Script summary

Phylogenetics benefits from using a large number of putatively independent nuclear loci and their combination with other sources of information, such as the plastid and mitochondrial genome. Selecting such orthologous low-copy nuclear (LCN) loci is still a challenge for non-model organisms. In recently published phylogenies based on target enrichment of several hundred LCN genes, these loci were selected from transcriptomes, genomes, gene expression studies, the literature, or a combination of these sources. Automated bioinformatic pipelines for the selection of LCN genes are, however, largely absent. We created a user-friendly, automated and interactive script named Sondovač to design LCN loci by a comparison between transcriptome and genome skim data. The script is licensed under open-source license GPL v.3 allowing further modifications. It runs on major Linux distributions and Mac OS X. Strong bioinformatics skills and access to high-performance computer clusters are not required; Sondovač runs on a standard desktop computer equipped with modern CPU.

# Pipeline --- how the data are processed

A transcriptome assembly and paired-end genome skim raw data are combined to get hundreds of orthologous LCN loci. Enrichment of multi-copy loci is minimized by using unique transcripts only, which are obtained by comparing all transcripts and removing those sharing ≥90% sequence similarity using BLAT. Before matching the genome skim data against those unique transcripts, reads of plastid (and mitochondrial) origin are removed with Bowtie 2 and SAMtools utilizing reference sequences. Paired-end reads are subsequently combined with FLASH. These processed reads are matched against the unique transcripts sharing ≥85% sequence similarity with BLAT. Transcripts with >1000 BLAT hits (indicating repetitive elements) and BLAT hits containing masked nucleotides are removed before de novo assembly of the BLAT hits to larger contigs with Geneious, using the medium sensitivity / fast setting. After assembly, only those contigs that comprise exons of a minimum bait length (usually ≥120 bp in case of probe design for phylogenies) and have a certain minimum total locus length (multiple of the bait length, should not be too short in order to obtain sufficient phylogenetically informative signal; we recommend at least ≥600 bp) are retained. To ensure that probes do not target multiple similar loci, any probe sequences sharing ≥90% sequence similarity are removed using cd-hit-est, followed by a second filtering step for contigs containing exons of a minimum bait length and totaling minimum loci length (see comments above). To ensure that plastid sequences are absent from the probes, the probe sequences are matched against the plastome reference sharing ≥90% sequence similarity with BLAT and the hits removed from the probe set. The workflow of Sondovač is summarized in Figure 1 of the PDF manual. Sondovač has three parts: two script parts and an intermediate part using Geneious. Input files for `sondovac_part_a.sh` are FASTA transcriptome data, FASTQ paired-end genome skim reads and a plastome (and possibly also mitochondriome) reference. The input file for Geneious is the output of `sondovac_part_a.sh`. The output files of Geneious are input files for `sondovac_part_b.sh`. The output file of `sondovac_part_b.sh` is the final list of probes, from which the user needs to remove the putative plastid sequences indicated in step 11 of the pipeline. In case of detection of remaining plastid sequences you will have to remove those plastid sequences manually from the output file of sondovac_ part_b.sh (see detailed documentation of `sondovac_part_b.sh` in the PDF manual).

# Software dependencies and installation

## Requirements to run Sondovač

Sondovač uses several scientific software packages and UNIX tools. See `INSTALL.md` and PDF manual for details.

## Geneious

Geneious is a DNA alignment, assembly, and analysis software and one of the most common software platforms used in genomics. It is utilized for de novo assembly in Sondovač. We plan to replace it by a free open-source command line tool in a future release of Sondovač. Visit https://www.geneious.com/ for download, purchase, installation and usage of Geneious. After the input data are processed (interactively or not) by `sondovac_part_a.sh`, the user must process its output manually by Geneious according to the instructions given in the PDF manual. The output of Geneious is then processed by `sondovac_part_b.sh`, which produces the final probe set. Geneious was tested with versions 6, 7, 8 and 9. See PDF manual for details.

# Command-line parameters to run Sondovač

## General parameters

Shared by `sondovac_part_a.sh` as well as `sondovac_part_b.sh`.

* `-h`, `-v` --- Print help message and exit.
* `-l` --- Display `LICENSE.md` for license information (this script is licensed under GNU GPL v.3, other software under variable licenses). Exit viewing by pressing the "Q" key.
* `-r` --- Display `README.md` (this file) for detailed usage instructions. Exit viewing by pressing the "Q" key. More information is available in the  PDF manual.
* `-p` --- Display `INSTALL.md` for detailed installation instructions. Exit viewing by pressing the "Q" key. More information is available in the PDF manual.
* `-e` --- Display detailed citation information and exit. See the PDF manual for more information.
* `-o` --- Set name of output files. Output files will start with that name. Do not use spaces or special characters --- some software can not handle them correctly. Default value (if the user does not provide another name) is "output". See below for the list of produced output files.

## Input files

Please, use file names without spaces and without special characters.

* `-f FILE` --- Transcriptome input file in FASTA format.
  * Used in `sondovac_part_a.sh`
* `-c FILE` --- Plastome reference sequence input file in FASTA format.
  * Used in `sondovac_part_a.sh`, `sondovac_part_b.sh`
  * Plastome reference sequences from taxa up to the same order of the studied plant group are suitable. See Shannon C K. Straub, Matthew Parks, Kevin Weitemier, Mark Fishbein, Richard C. Cronn and Aaron Liston; American Journal of Botany (2012) 99(2): 349-364, <https://bsapubs.onlinelibrary.wiley.com/journal/15372197>
* `-m FILE` --- Mitochondriome reference sequence input file in FASTA format (optional)
  * Used in `sondovac_part_a.sh`
  * This step is optional, as plant mitochondrial genomes have largely variable sizes and high rearrangement rates.
* `-t FILE` --- Paired-end genome skim input file in FASTQ format (first file).
  * Used in `sondovac_part_a.sh`
* `-q FILE` --- Paired-end genome skim input file in FASTQ format (second file).
  * Used in `sondovac_part_a.sh`
* `-x FILE` --- Input file in TSV format (output of Geneious assembly).
  * Used in `sondovac_part_b.sh`
* `-z FILE` --- Input file in FASTA format (output of Geneious assembly).
  * Used in `sondovac_part_b.sh`

## Optional parameters

See chapter "Pipeline" in the PDF manual for steps referred here. If those parameters are not provided, the default values are used.

* `-a ###` --- Maximum overlap length expected in approximately ≥90% of read pairs (parameter -M of FLASH, see its manual for details). FLASH can not combine paired-end reads that do not overlap by at least 10 bp (default minimum overlap length).
  * Step 4 of Sondovač, `sondovac_part_a.sh`.
  * DEFAULT: 65
  * OPTIONS: Integer ranging from 10 to 300
* `-y ##` --- Sequence similarity between unique transcripts and the filtered, combined genome skim reads (parameter -minIdentity of BLAT, see its manual for details). Filtering for orthologs, using sequence similarity as criterion.
  * Step 5 of Sondovač, `sondovac_part_a.sh`.
  * DEFAULT: 85 (highly recommended)
  * OPTIONS: Integer ranging from 70 to 100
* `-g` --- Choice of transcript or genome skim sequences for further processing. Depending on the phylogenetic depth that should be obtained, the probe sequences need to be designed from either the transcript or genome skim sequences, or it might not matter (if the taxa, from which the transcriptome and genome skim data were generated, are closely related).
  * Step 6.1 of Sondovač, `sondovac_part_a.sh`.
  * DEFAULT: no usage of -g (genome skim sequences)
  * OPTIONS: usage of -g (transcript sequences)
* `-s ####` --- Number of BLAT hits per transcript when matching unique transcripts and the filtered, combined genome skim reads. Transcripts with a high number of BLAT hits, indicating repetitive elements, need to be removed from the putative probe sequences.
  * Step 6.2 of Sondovač, `sondovac_part_a.sh`.
  * DEFAULT: 1000
  * OPTIONS: Integer ranging from 100 to 10000
* `-b ###` --- Minimum exon (bait) length. The minimum exon length should not fall below the bait length in order to account for specific binding between genomic libraries and baits during hybridization.
  * Steps 8 and 10 of Sondovač, `sondovac_part_b.sh`.
  * DEFAULT: 120 (preferred length for phylogeny).
  * OPTIONS: 80, 100, 120
* `-k ###` --- Minimum total locus length. When running the script in interactive mode, the user will be asked which value to use. A table summarizing the total number of LCN loci, which will be the result of the probe design for all minimum total locus lengths that the user can select (600 bp, 720 bp, 840 bp, 960 bp, 1080 bp, 1200 bp), will be displayed to facilitate this choice.
  * Steps 8 and 10 of Sondovač, `sondovac_part_b.sh`.
  * DEFAULT: 600
  * OPTIONS: 600, 720, 840, 960, 1080, 1200
* `-d 0.##` --- Sequence similarity between probe sequences (parameter -c of cd-hit-est, see its manual for details). Probes that target multiple similar loci need to be removed.
  * Step 9 of Sondovač, `sondovac_part_b.sh`.
  * DEFAULT: 0.9 (highly recommended).
  * OPTIONS: Decimal ranging from 0.85 to 0.95
* `-y ##` --- Sequence similarity between probe sequences and plastome reference (parameter -minIdentity of BLAT, see its manual for details). Some plastid reads might not have been removed in step2; they should be removed by this step.
  * Step 11 of Sondovač, `sondovac_part_b.sh`.
  * DEFAULT: 90 (highly recommended)
  * OPTIONS: Integer ranging from 70 to 100

## Examples

The basic and most simple usage (running in interactive mode):

    ./`sondovac_part_a.sh` -i

Specify some of the required input files, otherwise run interactively:

    ./`sondovac_part_a.sh` -i -f input.fa -t reads1.fastq -q reads2.fastq

Running in non-interactive, automated mode:

    ./`sondovac_part_a.sh` -n -f input.fa -c referencecp.fa -m referencemt.fa -t reads1.fastq -q reads2.fastq

Modify parameter "-a", otherwise run interactively:

    ./`sondovac_part_a.sh` -i -a 300

Run in non-interactive mode (parameter "-n") --- in such case the user must specify all required input files (parameters "-f", "-c", "-m", "-t" and "-q"). Moreover, parameter "-y" is modified:

    ./`sondovac_part_a.sh` -n -f input.fa -c referencecp.fa -m referencemt.fa -t reads1.fastq -q reads2.fastq -y 90

Modifying parameter "-s". Note that the interactive mode -i is implicit and does not need to be specified explicitly:

    ./`sondovac_part_a.sh` -s 950

# Input and output files

All names of input files and paths to them must be without spaces and without special characters (some software have difficulties handling these).

**Script `sondovac_part_a.sh` requires as input files:**

- **1)** Transcriptome input file in FASTA format. Note: For technical reasons, the labels of FASTA sequences must be unique numbers (no other characters). Sondovač will check the labels, and if they are not in an appropriate form, a copy of this input file with correct labels will be created.
- **2)** Plastome reference sequence input file in FASTA format.
- **3)** Paired-end genome skim input file in FASTQ format (two files --- forward and reverse reads).
- **4)** OPTIONAL: Mitochondriome reference sequence input file in FASTA format. This file is not required.

**Script `sondovac_part_a.sh` creates the following files:**

- **1)** `*_renamed.fasta` --- If needed, copy of the transcriptome input file with the changed labels of the FASTA sequences (unique numbers corresponding to the line numbers in the  original file) will be created. File 
- **2)** `*_old_and_new_names.tsv` then contains two columns: 1) the original sequence labels as in the user-provided transcriptome input file and 2) new sequence labels. This might be useful to trace back certain sequences/probes. 1)  `*_blat_unique_transcripts.psl` --- Output of BLAT (removal of transcripts sharing ≥90% sequence similarity).
- **3)** `*_unique_transcripts.fasta` --- Unique transcripts in FASTA format.
- **4)** `*_genome_skim_data_no_cp_reads*` --- Genome skim data without cpDNA reads.
- **5)** `*_genome_skim_data_no_cp_no_mt_reads*` --- Genome skim data without mtDNA reads --- only if mitochondriome reference sequence was used.
- **6)** `*_combined_reads_co_cp_no_mt_reads*` --- Combined paired-end genome skim reads.
- **7)** `*_blat_unique_transcripts_versus_genome_skim_data.pslx` --- Output of BLAT (matching of the unique transcripts and the filtered, combined genome skim reads sharing ≥85% sequence similarity).
- **8)** `*_blat_unique_transcripts_versus_genome_skim_data.fasta` --- Matching sequences in FASTA.
- **9)** `*_blat_unique_transcripts_versus_genome_skim_data-no_missing_fin.fsa` --- Final FASTA sequences for usage in Geneious.

Files 1--8 are not necessary for further processing by this pipeline, but may be useful for the user. The last file (9) is used as input file for Geneious in the next step.

[Geneious](https://www.geneious.com/) requires as input the last output file of `sondovac_part_a.sh` (file 9: `*_blat_unique_transcripts_versus_genome_skim_data-no_missing_fin.fsa`).

Output of Geneious are two exported files (see Geneious chapter in the PDF manual):

- **1)** Final assembled sequences exported as TSV.
- **2)** Final assembled sequences exported as FASTA.

Script `sondovac_part_b.sh` requires as input files:

- **1)** Plastome reference sequence input file in FASTA format.
- **2)** Assembled sequences exported from Geneious as TSV.
- **3)** Assembled sequences exported from Geneious as FASTA.

Script `sondovac_part_b.sh` creates the following files:

- **1)** `*_prelim_probe_seq.fasta` --- Preliminary probe sequences.
- **2)** `*_prelim_probe_seq_cluster_100.fasta` --- Unclustered exons and clustered exons with 100% sequence identity.
- **3)** `*_prelim_probe_seq_cluster_90.clstr` --- Unclustered exons and clustered exons with more than a certain sequence similarity (CLSTR file).
- **4)** `*_unique_prelim_probe_seq.fasta` --- Unclustered exons / exons with less than a certain sequence similarity.
- **5)** `*_similarity_test.fasta` --- Contigs that comprise exons ≥ bait length and have a certain total locus length.
- **6)** `*_target_enrichment_probe_sequences_with_pt.fasta` --- All probes in FASTA, with putative plastid sequences (if there were any BLAT hits, putative plastid sequences are listed in next file).
- **7)** `*_possible_cp_dna_gene_in_probe_set.pslx` --- In case of any BLAT hits, putative re- maining plastid probe sequences from *_target_enrichment_probe _sequences_with_pt.fasta are listed here. Not removing plastid genes will take lots of space on the Illumina lane for enriched plastid reads that should actually be available for enriched nuclear reads.
- **8)** `*_target_enrichment_probe_sequences.fasta` --- Final probes in FASTA.

An asterisk (*) denotes the beginning of the output files' names specified by the user with parameter `-o`. If user does not select a custom name, default value ("output") will be used.

By default, output files are created in the same directory from which Sondovač was launched. Output files can be saved in a custom directory by specifying an output directory together with parameter `-o`:

```shell
# Find current directory (e.g. /home/user):
pwd
# Launching Sondovač located in directory /home/user/sondovac and save output to e.g. desktop (/home/user/Desktop):
./sondovac/`sondovac_part_a.sh` -o Desktop/MyFile
# Sondovač will save software (if needed) in "bin" directory located in directory from which it was launched, see it:
ls bin/*
# Output files are in desired directory, see them e.g. by:
ls -lh Desktop/MyFile*
```

# Sample data

Together with the script, we provide the ZIP archive (~1.8 GB) containing example input files for running the script: Oxalis genome skimming data and transcriptome as well as the Ricinus cpDNA and mtDNA reference sequences. See https://github.com/V-Z/sondovac/wiki/Sample-data for download of sample data.

# License

The script is licensed under permissive open-source license allowing redistribution and modifications. Check file `LICENSE.md` or PDF manual for details and feel free to enhance the script.

# Citation

When using Sondovač please cite:

Roswitha Schmickl, Aaron Liston, Vojtěch Zeisek, Kenneth Oberlander, Kevin Weitemier, Shannon C.K. Straub, Richard C. Cronn, Léanne L. Dreyer and Jan Suda
Phylogenetic marker development for target enrichment from transcriptome and genome skim data: the pipeline and its application in southern African Oxalis (Oxalidaceae)
Molecular Ecology Resources (2016)
http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12487/abstract

# Other questions not covered here and reporting problems

If you have a question or you encounter a problem, please see https://github.com/V-Z/sondovac/issues and feel free to ask any question and/or express any wish. The authors will do their best to help you.

