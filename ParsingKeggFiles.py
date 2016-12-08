#This script parses a variety of files associated with Kegg annotation
#of protein fasta files. It can parse the Ghost Koala annotation
#and also a simple mapping of pathways to KO accessions


#this is written in Python 3 

HelpOutput = """
python ParsingKeggFiles.py -k file.txt -m file.txt -n file.txt

Required Parameters:
-k   <file>  This is a file from Ghost Koala annotation.
             This is the output from the detailed download (six columns)
-m   <file>  This is the file accessions to gene mapping
-n   <file>  This is the file of gene to KO mapping

Example: python .\ParsingKeggFiles.py -k .\GhostKoalaAnnotation\Detailed_annotation -k .\kegg_desc.txt -o .\KeggAnnotationTable.txt
"""


import os
import sys
import getopt
import string




class CreatorClass:

    def __init__(self):
        "nothing to put in really"
        self.KoalaAnnotationDirectory = "" #this directory holds a bunch of files for annotating genomes
        self.ProteinAccessionToKO = {} # key = protein accession, value = KO
        self.PathwayKOAssociationFile = "" #this is a small
        self.KOsInAPathway = {} # key = pathway, value = list of KOs
        self.PathwayDescriptions = {} #key = pathway, value = description



    def Main(self):
        #run through all the Koala Annotations
        KoalaFileList = os.listdir(self.KoalaAnnotationDirectory)
        for FileName in KoalaFileList:
            #freaking windows thumbs.db
            if not FileName[-3:] == "txt":
                continue
            FullPath = os.path.join(self.KoalaAnnotationDirectory, FileName)
            self.ParseGenomeAnnotationFromKoala(FullPath)
        self.ParsePathwayAssociation(self.PathwayKOAssociationFile)


    def ConvertListOfAccessionsToKOs(self, ListOfAccessions):
        #just auxiliary method
        ToReturn = []
        for Acc in ListOfAccessions:
            if not Acc in self.ProteinAccessionToKO:
                continue
            KO = self.ProteinAccessionToKO[Acc]
            ToReturn.append(KO)
        return ToReturn


    def ParsePathwayAssociation(self, Path):
        #here I have a simple two column file that tells me which pathway a given KO belongs to
        """
        ko:K00001   path:map00010   Glycolysis / Gluconeogenesis
        ko:K00001   path:map00071   Fatty acid degradation
        ko:K00001   path:map00350   Tyrosine metabolism
        ko:K00001   path:map00625   Chloroalkane and chloroalkene degradation
        """
        Handle = open(Path, 'r')
        Header = Handle.readline() #pop it off.
        for Line in Handle:
            (KO_precursor, Path_precursor, Path_Description) = Line.strip().split("\t")
            KO = KO_precursor.replace("ko:", "")
            Pathway = Path_precursor.replace("path:", "")
            #now assign them to a pathway
            if not Pathway in self.KOsInAPathway:
                self.KOsInAPathway[Pathway] = [] #empty list
                self.PathwayDescriptions[Pathway] = Path_Description
            self.KOsInAPathway[Pathway].append(KO)

        Handle.close()




    def ParseGenomeAnnotationFromKoala(self, Path):
        #trying to get out the gene annotations of Kegg
        
        Handle= open(Path, 'r')
        for line in Handle: #python 3 a handle is natively iterable

            Bits = line.split("\t") #don't strip all the white space, because empty cells makes a \t\t\t string at the end which is important
            # [0] is identifier
            #[1] is the Kegg_ortholog if it passes the best cutoff
            #[2] is name of the kegg ortholog
            #[3] is a score
            #[4] is a Kegg_ortholog assignment that was lower quality
            #[5] is a score for the lower quality assignment, or something, I'm not totally sure
            KO = Bits[1]
            FastaLineIdentifier = Bits[0]  #sp|E1WTS4|BIOAB_BACF6
            #Accession = FastaLineIdentifier.split("|")[1]
            # I believe that what is parsed out by MSGF and Koala will be the same, so let's just leave it alone
            Accession = FastaLineIdentifier

            #check to see if this line has an annotation
            if KO == "":
                #we don't really care about proteins that don't get assigned to a KO. Let's leave them out
                continue
            #print ("Putting %s and %s into the dictionary"%(KO, Accession))
            self.ProteinAccessionToKO[Accession] = KO


        Handle.close()

       
    def AddFiles(self, _KoalaDir, _PathwayKOAssociationFile):
        #this is called by the iPython Notebook to mimic the ParseCommandLine() function
        #but from the notebook
        self.KoalaAnnotationDirectory = _KoalaDir
        self.PathwayKOAssociationFile = _PathwayKOAssociationFile


