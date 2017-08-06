#!/usr/bin/python

#Copyright 2017, Justin R. Klesmith
#All rights reserved.
"""Pipeline for building PSSMs for protein analysis"""

#Note: If you are defining custom regions to split the MSA for psiblast go to line 398

from __future__ import division
from argparse import ArgumentParser
from xml.etree import cElementTree
import os
import time

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2017, Justin R. Klesmith"
__license__ = "BSD-3"
__version__ = "2.1, Build: 20170806"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

#Setup subparsers
parser = ArgumentParser(prog='Protein PSSM Pipeline', description='Create a PSSM for a protein using a BlastP search')
subparsers = parser.add_subparsers(help='Usage: python ProteinPSSMPipeline.py {RunMode} -flags', dest='RunUsed')

# Create the parser for the XMLtoFASTA run
parser_a = subparsers.add_parser('XMLtoFASTA', help='Convert the XML file into a fasta file for CDHIT')
parser_a.add_argument('-f', dest='ncbi', action='store', required=True, help='File path to XML file from the BlastP search')
parser_a.add_argument('-w', dest='wildtype', action='store', required=True, help='Path to the file with the wild-type sequence')
parser_a.add_argument('-o', dest='outfile', action='store', required=True, help='Output filename')
parser_a.add_argument('-q', dest='querylen', action='store', nargs='?', const=1, default=0.6, help='Minimum length of match to query length. Default = 0.6 (i.e. 60 percent)')
parser_a.add_argument('-s', dest='seqid', action='store', nargs='?', const=1, default=0.3, help='Minimum sequence identity of hit. Default = 0.3 (i.e. 30 percent)')

# Create the parser for the CDHITWildTypeCheck run
parser_b = subparsers.add_parser('CDHITWTCheck', help='Checks to see if WT is still there after CDHIT')
parser_b.add_argument('-a', dest='afafile', action='store', required=True, help='File path to the CDHIT output file')
parser_b.add_argument('-w', dest='wildtype', action='store', required=True, help='Path to the file with the wild-type sequence')

# Create the parser for the ProcessMSA run
parser_c = subparsers.add_parser('ProcessMSA', help='Convert the MSA for use with PSI-Blast')
parser_c.add_argument('-m', dest='msa', action='store', required=True, help='MSA file path')
parser_c.add_argument('-w', dest='wildtype', action='store', required=True, help='Path to the file with the wild-type sequence')
parser_c.add_argument('-n', dest='maxhits', action='store', nargs='?', help='Maximum number of hits to send to PSI-Blast (typically 500)')

# Create the parser for the SplitAlignmentForPSIBlast run
parser_d = subparsers.add_parser('SplitForPSIBlast', help='Split the processed MSA')
parser_d.add_argument('-m', dest='msa', action='store', required=True, help='MSA file path')
parser_d.add_argument('-n', dest='filename', action='store', required=True, help='Name of output files')
parser_d.add_argument('-s', dest='SizeOfRegions', action='store', nargs='?', help='Size of the regions')
parser_d.add_argument('-w', dest='wildtype', action='store', required=True, help='Path to the file with the wild-type sequence')

# Create the parser for the PSIBlast command line prep run
parser_e = subparsers.add_parser('PSIBlastCmdLine', help='Prepare a shell script file for PSIBlast')
parser_e.add_argument('-p', dest='path', action='store', required=True, help='PSIBlast file path')
parser_e.add_argument('-n', dest='name', action='store', required=True, help='Name of subalignments ie name given to the splitalignment step')
parser_e.add_argument('-w', dest='wildtype', action='store', required=True, help='Path to the file with the wild-type sequence')

# Create the parser for combining the PSSM outputs into a CSV heatmap
parser_f = subparsers.add_parser('PSSMHeatmap', help='Combine the PSSM outputs and then make a CSV heatmap')
parser_f.add_argument('-o', dest='outname', action='store', required=True, help='Output filename')
parser_f.add_argument('-w', dest='wildtype', action='store', required=True, help='Path to the file with the wild-type sequence')

# parse some argument lists
args = parser.parse_args()

def WildType():
    """Get the wild-type protein information from the given fasta file"""
    f = open(args.wildtype, 'r')
    lines = f.readlines()
    WTNAME = lines[0].rstrip('\r\n')
    WTSEQ = lines[1].rstrip('\r\n')
    WTLEN = len(WTSEQ)
    f.close()
    return WTNAME, WTSEQ, WTLEN

def CleanStrForFilenames(filename):
    """Sanitize a string to be used as a filename"""
    keepcharacters = (' ', '.', '_')
    FILEOUT = "".join(c for c in filename if c.isalnum() or c in keepcharacters).rstrip()
    return FILEOUT

def XMLtoFASTA():
    """Convert the XML file into a fasta file for CDHIT"""
    #Check inputs for valid values
    try:
        MINLEN = float(args.querylen)
    except ValueError:
        print "XMLtoFASTA Error: An incorrect value given for the querylen (option -q)"
        quit()

    try:
        MINSEQID = float(args.seqid)
    except ValueError:
        print "XMLtoFASTA Error: An incorrect value given for the seqid (option -s)"
        quit()

    if MINLEN > 1:
        print "XMLtoFASTA Error: The minimum lenth to query length is greater than one. (Valid values: 0.0 to 1.0)"
        quit()
    elif MINLEN < 0:
        print "XMLtoFASTA Error: The minimum length to query length is less than zero. (Valid values: 0.0 to 1.0)"
        quit()

    if MINSEQID > 1:
        print "XMLtoFASTA Error: The minimum sequence identity is greater than one. (Valid values: 0.0 to 1.0)"
        quit()
    elif MINSEQID < 0:
        print "XMLtoFASTA Error: The minimum sequence identity is less than zero. (Valid values: 0.0 to 1.0)"
        quit()

    if os.path.isfile(args.wildtype) == False:
        print "XMLtoFASTA Error: The path to the file with the wild-type sequence is missing"
        quit()

    if os.path.isfile(args.ncbi) == False:
        print "XMLtoFASTA Error: Cannot open the processed NCBI TSV"
        quit()

    #Get the Wild-Type amino acid sequence
    WTNAME, WTSEQ, WTLEN = WildType()
    
    #Check to see if the output filename is valid
    OUTFILENAME = CleanStrForFilenames(args.outfile)
    
    #Flag to see if our WT is in the XML input
    WTINNCBI = False 

    #Import NCBI information from TSV
    f2 = open(OUTFILENAME, 'w')
    for event, elem in cElementTree.iterparse(args.ncbi):
        if elem.tag == "Hit":
            hitNUM = elem.find("Hit_num").text
            hitID = elem.find("Hit_id").text
            hitNAME = elem.find("Hit_def").text
            hitACCESSION = elem.find("Hit_accession").text
            hitIDENTITIES = int(elem.find("Hit_hsps").find("Hsp").find("Hsp_identity").text)
            hitALIGNLEN = int(elem.find("Hit_hsps").find("Hsp").find("Hsp_align-len").text)
            hitSEQUENCE = elem.find("Hit_hsps").find("Hsp").find("Hsp_hseq").text
            
            #Verify that the alignment length is >= 60% of query
            if float(hitALIGNLEN/WTLEN) >= MINLEN:

                #Verify that the sequence identity is >= 30%
                if float(hitIDENTITIES/hitALIGNLEN) >= MINSEQID:

                    #Remove dashes from the ncbi match (cd-hit shits them out as bad formatting)
                    chars_to_remove = ['-']

                    #Check to see if the sequence was WT
                    if hitSEQUENCE.translate(None, ''.join(chars_to_remove)) == WTSEQ:
                        WTINNCBI = True

                        #Write the wild-type to the file
                        f2.write(WTNAME+'\n')
                        f2.write(hitSEQUENCE.translate(None, ''.join(chars_to_remove))+'\n\n')
                    else:
                        #Write the other matches to the file
                        f2.write(">"+hitACCESSION+"_"+hitID+"_"+hitNAME+'\n')
                        f2.write(hitSEQUENCE.translate(None, ''.join(chars_to_remove))+'\n\n')
            
            elem.clear()

    #If WT was not in the NCBI input then add it
    if WTINNCBI == False:
        f2.write(WTNAME+'\n')
        f2.write(WTSEQ+'\n')

    #Close the output file
    f2.close()

    return

def CDHITWTCheck():
    """Checks to see if WT is still there after CDHIT"""
    #Check inputs for valid values
    if os.path.isfile(args.wildtype) == False:
        print "CDHITWTCheck Error: The path to the file with the wild-type sequence is missing"
        quit()

    if os.path.isfile(args.afafile) == False:
        print "CDHITWTCheck Error: Cannot open the CDHIT output file"
        quit()

    #Get the Wild-Type amino acid sequence
    WTNAME, WTSEQ, WTLEN = WildType()
    
    #Flag to see if our WT is in the input
    WTINNCBI = False

    #Open the file and check it line by line
    with open(args.afafile, 'r') as infile: #Open the file with the wild-type protein sequence
        for line in infile:

            #Check to see if the sequence was WT
            if line.rstrip('\n\r') == WTNAME:
                WTINNCBI = True

    #If WT was not in the NCBI input then add it
    if WTINNCBI == False:
        f2 = open(args.afafile, 'a')
        f2.write(WTNAME+'\r')
        f2.write('\n')
        f2.write(WTSEQ+'\r')
        f2.close()

    return

def ProcessMSA():
    """Convert the MSA for use with PSI-Blast"""
    #Check the inputs
    if os.path.isfile(args.wildtype) == False:
        print "ProcessMSA Error: The path to the file with the wild-type sequence is missing"
        quit()

    if os.path.isfile(args.msa) == False:
        print "ProcessMSA Error: Cannot open the msa file"
        quit()

    if args.maxhits == None:
        print "ProcessMSA: The maximum number of hits will be uncapped"
        MAXHITS = -1
    else:
        try:
            MAXHITS = int(args.maxhits)
            if MAXHITS < 1:
                print "ProcessMSA Error: The maximum number of hits will be uncapped due to a negative or zero value given"
                MAXHITS = -1
        except ValueError:
            print "ProcessMSA Error: The maximum number of hits will be uncapped due to an error on the command line argument given"
            MAXHITS = -1

    #Get the Wild-Type amino acid sequence
    WTNAME, WTSEQ, WTLEN = WildType()

    #
    #Step one: Import msa alignment from MUSCLE and make it one line per sequence
    Alignment = ""
    Output = ""
    with open(args.msa, 'r') as infile: #Open the file with the wild-type protein sequence   
        for line in infile:
            #Check to see if we have a header
            if line[0] == ">":

                if len(Output) > 0: #Ignores empty output
                    Alignment = Alignment + Output + "\n" #Add the current output to the growing alignment varible

                Output = "" #Empty the current alignment
                Output = Output + line.rstrip('\n') + "," #Assemble the first line of the new sequence
            else:
                Output = Output + line.rstrip('\n') #Keep assembling the line
    f = open('msatemp.csv', 'w')
    f.write(Alignment)
    f.close()

    #
    #Step two: Import MSA into a lookup table
    MSATable = {}
    Output = ""
    with open('msatemp.csv', 'r') as infile: #Open the previous file
        for line in infile:
            NUMCOMMAS = line.count(',') #Count the number of commas (the fix for multiple commas in name)
            split = line.split(",")
            if len(line) > 10: #Avoids empty lines from the file
                MSATable.update({split[0] : split[NUMCOMMAS].rstrip("\n")}) #Add a new entry to the array
                Output = Output + split[NUMCOMMAS] #Write the sequence (the highest numbered split entry)
    f = open('msatemp2.csv', 'w')
    f.write(Output)
    f.close()

    #
    #Step three: Mark the insertions with the letter Z (fixed from X as X is also used in the blast hits)
    ZedOut = ""
    WildtypeMSA = MSATable[WTNAME] + "\n" #Get the wild-type MSA sequence (with insertions)
    WTMSALEN = len(WildtypeMSA)
    with open('msatemp2.csv', 'r') as infile:
        for line in infile:
            for i in xrange(0, WTMSALEN):
                if WildtypeMSA[i] == "-":
                    ZedOut = ZedOut + "Z"
                else:
                    ZedOut = ZedOut + line[i]
    f = open('msatemp3.csv', 'w')
    f.write(ZedOut)
    f.close()

    #
    #Step four: Delete the insertions
    Output = ""
    with open('msatemp3.csv', 'r') as infile:
        for line in infile:
            Len = len(line)
            for i in xrange(0, Len):
                if line[i] != "Z":
                    Output = Output + line[i]
    f = open('msatemp4.csv', 'w')
    f.write(Output)
    f.close()

    #
    #Step five: Put wild-type on top, re-order the sequences by completeness, and cap the hits for psi-blast
    PSITable = []
    with open('msatemp4.csv', 'r') as infile:
        for line in infile:

            #Set the WT to the top and then add the rest into a list with their counts of dashes
            if line == WTSEQ + "\n":
                PSITable.insert(0, {'sequence' : line, 'counts' : line.count('-'), 'wt' : True}) #Move the wild-type to the top
            else:
                PSITable.append({'sequence' : line, 'counts' : line.count('-'), 'wt' : False})

    #Sort the table into a new list by counts and wt
    PSITable.sort(key=lambda k : (k['counts'], -k['wt']))

    #Iterate the list to a writeable string
    Output = ""
    Count = 1
    for x in PSITable:
        if MAXHITS == -1: #Unlimited number of output sequences
            Output = Output + x['sequence']
        else:
            if Count <= MAXHITS: #Add to the output if below or equal to the max given
                Output = Output + x['sequence']
                Count = Count + 1

    f = open('MSAForSplitting.csv', 'w')
    f.write(Output)
    f.close()

    return
    
def SplitForPSIBlast():
    """Split the processed MSA"""
    #Check the inputs
    if os.path.isfile(args.wildtype) == False:
        print "SplitForPSIBlast Error: The path to the file with the wild-type sequence is missing"
        quit()

    if os.path.isfile(args.msa) == False:
        print "SplitForPSIBlast Error: Cannot open the msa file"
        quit()

    if args.SizeOfRegions == None:
        print "SplitForPSIBlast: The number of regions is not defined, probably set manually?"
        RegionSize = None
    else:
        try:
            RegionSize = int(args.SizeOfRegions)
        except:
            print "SplitForPSIBlast Error: An incorrect input given for the size of regions"
            quit()
    if RegionSize < 1:
        print "SplitForPSIBlast Error: An incorrect input given for the size of regions"
        quit()
        
    #Sanitize the output filename
    FILEOUT = CleanStrForFilenames(args.filename)

    #Get the Wild-Type amino acid sequence
    WTNAME, WTSEQ, WTLEN = WildType()
    
    Regions = None

    #Create the regions from the given inputs
    if RegionSize >= 1:
        Regions = []
        RegionMod = divmod(WTLEN, RegionSize)
        WholeRegions = RegionMod[0]
        Remainder = RegionMod[1]

        #Add one to the total number of regions if there is a remainder
        if Remainder > 0:
            NumRegions = WholeRegions + 1
        else:
            NumRegions = WholeRegions

        #Calculate the low and high bounds of the regions
        for x in xrange(0, WholeRegions + 1):

            if x < WholeRegions:
                low = x * RegionSize
                high = (x * RegionSize) + RegionSize - 1
            else:
                low = x * RegionSize
                high = WTLEN - 1

            #Create the list from the regions
            Regions.append([low, high])

    print "SplitForPSIBlast: Calculated Regions:"
    print Regions
    print "SplitForPSIBlast: Total number of regions = " + str(NumRegions)

    #Assign sequence regions manually (this is after the X marked insertions are deleted)
    #Regions = [[[[0, 19], [20, 39], [40, 59], [60, 79], [80, 99], [100, 119], [120, 139], [140, 159], [160, 179], [180, 199], [200, 219], [220, 239], [240, 259], [260, 271]]]
    #NumRegions = 20

    #Final Check to see if our varibles are defined
    if Regions == None:
        print "SplitForPSIBlast Error: Regions is not defined"
        quit()

    if NumRegions == None:
        print "SplitForPSIBlast Error: The number of regions is not defined"
        quit()

    #Lets open the MSA (one sequence per line)
    with open(args.msa, 'r') as infile: #Open the file with the wild-type protein sequence
        for region in Regions:
            Start = region[0]
            End = region[1]
            counter = 0
            infile.seek(0)
            outfile = open('subalignment_'+FILEOUT+'_AAstart_'+str(Start)+'.fasta', 'w')
            for line in infile:
                Subline = ""
                for i in xrange(Start, End+1):
                    Subline = Subline + line[i].rstrip('\n')

                #Remove sequences with insertions
                if Subline.find('-') == -1:
                    outfile.write(">" + str(counter) + "\n")
                    outfile.write(Subline + "\n")
                    counter = counter + 1

            outfile.close()

    return
    
def PSIBlastCmdLine():
    """Prepare a shell script file for PSIBlast"""
    #Check the inputs
    if os.path.isfile(args.wildtype) == False:
        print "PSIBlastCmdLine Error: The path to the file with the wild-type sequence is missing"
        quit()

    if os.path.isfile(args.path) == False:
        print "PSIBlastCmdLine Error: Cannot open the msa file"
        quit()
        
    #Open the output file and then make the commands
    f = open('PSSMCommands.sh', 'w')
    for filename in os.listdir("."):
        if filename.startswith("subalignment_"+args.name+"_AAstart_"):
            file2 = filename.split('_')
            file3 = file2[3].split('.')
            SplitIndex = file3[0]
            f.write(args.path+" -subject "+args.wildtype+" -in_msa "+filename+" -out_ascii_pssm "+args.name+"_PSSM_"+SplitIndex+".pssm > "+args.name+"_PSSM_"+SplitIndex+".log\n")  
    f.close()
    return
    
def PSSMHeatmap():
    """Combine the PSSM outputs and then make a CSV heatmap"""

    #Sanitize the output filename
    FILEOUT = CleanStrForFilenames(args.outname)

    #Get the Wild-Type amino acid sequence
    WTNAME, WTSEQ, WTLEN = WildType()
    
    #The standard heatmap AA mappings
    AA_Table = 'FWYPMILVAGCSTNQDEHKR'
    
    #Based on the PSSM output file. First two columns are wt and location. The rest is the other columns with 0 = pssm score, 1 = percents
    PSSMTable = [['NoneZ', None], ['NoneZ', None], 
    ['A', 0], ['R', 0], ['N', 0], ['D', 0], ['C', 0], ['Q', 0], ['E', 0], ['G', 0], ['H', 0], ['I', 0], 
    ['L', 0], ['K', 0], ['M', 0], ['F', 0], ['P', 0], ['S', 0], ['T', 0], ['W', 0], ['Y', 0], ['V', 0], 
    ['A', 1], ['R', 1], ['N', 1], ['D', 1], ['C', 1], ['Q', 1], ['E', 1], ['G', 1], ['H', 1], ['I', 1], 
    ['L', 1], ['K', 1], ['M', 1], ['F', 1], ['P', 1], ['S', 1], ['T', 1], ['W', 1], ['Y', 1], ['V', 1]]
    Mutations = {} #Mutations matrix

    #Populate mutation matrix with None data
    for j in xrange(1, WTLEN+1):
        for i in enumerate(AA_Table):
            try:
                #Mutations[ResID][MutID[1]][0 = PSSMMetric, 1 = Percents]
                Mutations[j][i[1]] = [None, None]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None]

    #Read all of the files
    for filename in os.listdir("."):
        if filename.endswith(".pssm"):
            #Get the total number of lines in the PSSM file
            TOTALLINES = open(filename, 'r').read().count("\n") + 1

            #Get the starting index number
            file2 = filename.split('_')
            file3 = file2[2].split('.')
            SplitIndex = int(file3[0])

            #Open the file again for reading
            f = open(filename, 'r')
            lines = f.readlines()
            for i in xrange(3, TOTALLINES - 5):
                linesplit = lines[i].split(' ')
                ls2 = filter(None, linesplit)
                for j in xrange(2, 42):
                    Mutations[int(ls2[0])+SplitIndex][PSSMTable[j][0]][PSSMTable[j][1]] = ls2[j]
            f.close()

    #This makes a CSV style report of rows of letters and columns of residues
    #Print off the Number
    Numbering = " "
    for q in xrange(1, WTLEN+1):
        Numbering = Numbering+","+str(q)

    #Print off the WT Residue
    WTResi = " "
    for w in xrange(0, WTLEN):
        WTResi = WTResi+","+WTSEQ[w]

    #Print off the mutations
    Output = "Last position-specific scoring matrix computed\n"
    for i in enumerate(AA_Table):
        Output = Output+i[1]+","
        for j in xrange(1, WTLEN+1):
            Output = Output+str(Mutations[j][i[1]][0])+","
        Output = Output+"\n"

    Output = Output+"\nWeighted observed percentages rounded down\n"

    for i in enumerate(AA_Table):
        Output = Output+i[1]+","
        for j in xrange(1, WTLEN+1):
            Output = Output+str(Mutations[j][i[1]][1])+","
        Output = Output+"\n"

    #Write the heatmap to a newfile
    outfile = open('PSSM_Heatmap_'+FILEOUT+'.csv', 'w')
    outfile.write(Numbering+'\n')
    outfile.write(WTResi+'\n')
    outfile.write(Output)
    outfile.close()

    return

def main():
    """The main program subroutine"""
    print "ProteinPSSMPipeline"
    print "Author: "+__author__
    print "Contact: "+__email__[0]+", "+__email__[1]
    print __copyright__
    print "Version: "+__version__
    print "License: "+__license__
    print ""
    print "Please cite:"
    print "Github [user: JKlesmith] (www.github.com)"
    print "Klesmith, J. R., Bacik, J.-P., Wrenbeck, E. E., Michalczyk, R., and Whitehead, T. A. (2017) Trade-offs between enzyme fitness and solubility illuminated by deep mutational scanning, Proceedings of the National Academy of Sciences 114, 2265-2270."
    print ""
    print "Run parameters:"
    print time.strftime("%H:%M:%S")
    print time.strftime("%m/%d/%Y")
    print args
    
    #Select which run off of the command line args
    if args.RunUsed == "XMLtoFASTA":
        print "Run Mode: XMLtoFASTA"
        XMLtoFASTA()
    elif args.RunUsed == "CDHITWTCheck":
        print "Run Mode: CDHITWTCheck"
        CDHITWTCheck()
    elif args.RunUsed == "ProcessMSA":
        print "Run Mode: ProcessMSA"
        ProcessMSA()
    elif args.RunUsed == "SplitForPSIBlast":
        print "Run Mode: PSplitForPSIBlast"
        SplitForPSIBlast()
    elif args.RunUsed == "PSIBlastCmdLine":
        print "Run Mode: PSIBlastCmdLine"
        PSIBlastCmdLine()
    elif args.RunUsed == "PSSMHeatmap":
        print "Run Mode: PSSMHeatmap"
        PSSMHeatmap()

if __name__ == '__main__':
    main()
