#!/usr/bin/env python
from os import path
from optparse import OptionParser # Import the option parser module

###### OPTIONS and USAGE ######
parser = OptionParser(usage = """grab_singleton_clusters.py -i INFILE -o OUTFILE
grab_singleton_clusters.py -- Finds those clusters from a CD-HIT (and related
    programs) .clstr file that include only one sequence or multiple sequences
    with 100% identity (in which case it chooses the longest sequence). It
    produces a file in .clstr format. 
Copyright (c) 2014 Kevin Weitemier
Version 1.00
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. A copy of this license is available at <http://www.gnu.org/licenses/>.
Great effort has been taken to make this software perform its said
task, however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
If you find this program useful, please cite:
K. Weitemier, S.C.K. Straub, R. Cronn, M. Fishbein, R. Schmickl, A. McDonnell,
and A. Liston. 2014. Hyb-Seq: Combining target enrichment and genome skimming
for plant phylogenomics. Applications in Plant Sciences 2(9): 1400042.
Input - A .clstr file created by CD-HIT or related programs.""")
parser.add_option("-i", action="store", type="string", dest="inname",
    help="Input filename", default="")
parser.add_option("-o", action="store", type="string", dest="outname",
    help="Output filename", default="")
(options, args) = parser.parse_args()

# Make sure all filenames are given
if options.inname == "":
    parser.error("Please include an input file using -i.")
if options.outname == "":
    parser.error("Please include an output file using -o.")
if path.isfile(options.outname):
    parser.error("""The output filename exists. Please delete it first or \
choose a new name.""")

###### OPENING INPUT/OUTPUT FILES ######
Infile = open(options.inname, 'r')
Outfile = open(options.outname, 'w')

###### MAIN PROGRAM ######
Line = Infile.readline()

while Line:
    if Line.startswith('>'):
        StoredCluster = Line.strip()
        Line = Infile.readline()
        FirstName = Line.strip()
        Line = Infile.readline()
        if Line.startswith('>'):
            Outfile.write("%s\n%s\n" % (StoredCluster,FirstName))
#            continue
        else:
            if not Line:
                break
            else:
                ListOfLines, ListOfPercents = ([], [])
                ListOfLines.append(FirstName)
                while not Line.startswith('>'):
                    ListOfLines.append(Line.strip())
                    Line = Infile.readline()
                    if not Line:
                        break
                for Name in ListOfLines:
                    if 'at' in Name:
                        Fields = Name.split('/')
                        Percent = Fields[1].strip('%')
                        if Percent != '100.00':
                            break
                else:
                    LineToPrint = ""
                    for Alternate in ListOfLines:
                        if Alternate.endswith("*"):
                            LineToPrint = Alternate
                            break
                    Outfile.write("%s\n%s\n" % (StoredCluster,LineToPrint))

    else:
        Outfile.write("Got mis-placed, moving to next cluster!\n%s\n" % (Line))
        while not Line.startswith('>'):
            Line = Infile.readline()
#        continue

Outfile.write("%s\n%s\n" % (StoredCluster,FirstName))
###### CLEANUP ######
Infile.close()
Outfile.close()
# EOF
