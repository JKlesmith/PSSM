#!/bin/sh

#PSSM Pipeline
#Copyright 2017, Justin R. Klesmith
#BSD-3 License
#Github [user: JKlesmith] (www.github.com)
#Original procedures are outlined in:
#Klesmith, J. R., Bacik, J.-P., Wrenbeck, E. E., Michalczyk, R., and Whitehead, T. A. (2017) Trade-offs between enzyme fitness and solubility illuminated by deep mutational scanning, Proceedings of the National Academy of Sciences 114, 2265-2270.
#Version 2 Pipeline is outlined in:
#Publication in preparation

#You will need to make a file with the wild-type name and sequence (example below and in the archive)
#>WildType_LGK
#MPIATSTGDNVLDFTVLG...

#Convert the blastp xml output to a tsv file
#Edit: Add the XML filename
python BlastXMLtoTSV.py -c ext BLASTPsearchresult.xml > BlastP.tsv

#Convert the TSV into a fasta file with filters (-q is the min length of match sequence to query sequence, -s is the min sequence identity, -o is output filename for CDHIT)
#Edit: Add the path to the file with the wild-type sequence, and other optional paramaters
python ProteinPSSMPipeline.py TSVtoFASTA -f BlastP.tsv -q 0.6 -s 0.3 -w WildType -o Protein.fa

#Run CDHIT (-c sets how similar to cluster the higher the number = the higher the threshold to cluster, -T is the number of threads)
#Edit: Add the correct path to CD-HIT
.\Programs\Windows\cd-hit-v4.6.7-2017-0501\cd-hit.exe -i Protein.fa -o Protein.afa -c 0.995 -M 40000 -T 2

#WT Check - Checks to see if the wild-type sequence is still there after CD-HIT, adds it if not found
#Edit: Add the path to the file with the wild-type sequence
python ProteinPSSMPipeline.py CDHITWTCheck -a Protein.afa -w Wildtype

#Run MUSCLE
#Edit: Add the correct path to MUSCLE
.\Programs\Windows\muscle3.8.31_i86win32.exe -in Protein.afa -out Protein.msa

#Remove insertions from the MSA (-m is the msa filename from muscle, -n is the maximum number of seuqnces to send to psi-blast (for unlimited don't add the -n flag))
#Edit: Add the path to the file with the wild-type sequence
python ProteinPSSMPipeline.py ProcessMSA -m Protein.msa -w WildType -n 500

#Split the processed MSA
#Edit: Add the size of each region and path to the wildtype sequence (Note:custom regions can be set in the script manually), and the project name
python ProteinPSSMPipeline.py SplitForPSIBlast -m MSAForSplitting.csv -s SizeOfEachRegion -w WildType -n ProjectName

#Create the command lines for PSI-Blast
#Edit: Add the path to PSIBlast Program and the path to a file containing the wild-type sequence on a single line, and the project name
python ProteinPSSMPipeline.py PSIBlastCmdLine -p PathToPSIBlast -w WildType -n ProjectName

#Run PSI-Blast
sh PSSMCommands.sh

#Convert the PSSM outputs into a heatmap and then csv
#Edit: Add the name of the output file
python ProteinPSSMPipeline.py PSSMHeatmap -o OutfileName -w WildType
