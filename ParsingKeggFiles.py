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



class Protein:
    def __init__(self, accession):
        self.accession = accession
        self.ManualKO = ""
        self.KoalaKO = ""
    def addManualKO(self, KO):
        self.ManualKO = KO
    def addKoalaKO(self, KO):
        self.KoalaKO = KO
    def GetLine(self):
        Line = "%s\t%s\t%s"%(self.accession, self.KoalaKO, self.ManualKO)
        return Line


class CreatorClass:

    def __init__(self):
        "nothing to put in really"
        self.KoalaAnnotationPath = "" # this file holds GhostKoala kegg annotations. 
        self.ManualGeneKOPath = ""
        self.ManualGeneAccessionPath = ""
        self.GeneMappingDictionary = {} #key =gene, value = accession 
        self.IncludeAllAssignments = 0 #kegg has assignments that don't pass the top cutoff. so this is a flag
        #								# to include or exclude those
        self.ProteinObjectsDictionary = {} #key = Accession, value= object
        self.OutputPath = "renameYourOutput.txt"

    def Main(self):

        self.ParseGenomeAnnotationFromKoala(self.KoalaAnnotationPath)
        self.ParseGeneToAccession(self.ManualGeneAccessionPath)
        self.ParseKO_to_Gene(self.ManualGeneKOPath)
        #now I would need to loop through this and find those passing my cutoff
        self.printItAllOff()

    def ParseKeggFile(self, Path):
        #here our goal is to create a dictionary so that we can print off names when we say the KO
        Handle = open(Path, 'r')
        Header = Handle.readline() # pop it off. I don't really need it. but don't want it in the loop
        for line in Handle:
            #bits[0] is the KO
            #bits[1] is the gene name (unnessary)
            #bits[2] is the EC number
            #bits[3] is the definition, or plain english description of the ortholog family function
            Bits =line.split("\t")
            KO_identifier =Bits[0]
            KO_description = Bits[3].strip() # to get rid of the return line
            self.KeggOrthologDescriptions[KO_identifier] = KO_description
        Handle.close()

    def ParseKO_to_Gene(self, Path):
        """ here i'm trying to get out the KO assigned by manual curation and put
        that in the right object so that it is able to compare to Koala
        bfg:BF638R_0816 ko:K05349
        bfg:BF638R_0681 ko:K17103
        bfg:BF638R_3409 ko:K01843

        """
        Handle = open(Path, 'r')
        for Line in Handle:
            #[0] is the Gene
            #[1] is the KeggOrtholog identifier
            (Gene, KO) = Line.strip().split("\t")
            Gene = Gene.replace("bfg:", "")
            KO = KO.replace("ko:", "")
            if not Gene in self.GeneMappingDictionary:
                continue #skip this round, because it didn't make it into having a protein accession, probably just a tRNA
            Accession = self.GeneMappingDictionary[Gene]
            if Accession in self.ProteinObjectsDictionary:
                self.ProteinObjectsDictionary[Accession].addManualKO(KO)
            else:
                NewObject = Protein(Accession)
                NewObject.addManualKO(KO)
                self.ProteinObjectsDictionary[Accession] = NewObject
    

    def ParseGeneToAccession(self, Path):
        """Here I'm trying to store for a moment these types of relationships
        bfg:BF638R_0001 up:E1WJP6
        bfg:BF638R_0002 up:E1WJP7
        bfg:BF638R_0003 up:E1WJP8

        """

        Handle = open(Path, 'r')
        for Line in Handle:
            #[0] is the gene
            #[1] is the accession
            (Gene, Accession) = Line.strip().split("\t")
            #now strip out the silly  bfg: and up:
            Gene = Gene.replace("bfg:", "")
            Accession = Accession.replace("up:", "")
            self.GeneMappingDictionary[Gene] = Accession
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
            Accession = FastaLineIdentifier.split("|")[1]

            #check to see if this line has an annotation
            if KO == "":
                #now if directed try the secondary annotation
                if not self.IncludeAllAssignments:
                    continue
                KO = Bits[4]
                if KO == "":
                    continue
                #the implicit else means: we are allowed to include alternate assignments AND there is an alternate assignment

            NewObject = Protein(Accession)
            NewObject.addKoalaKO(KO)
            self.ProteinObjectsDictionary[Accession] = NewObject

        Handle.close()

    def printItAllOff(self):
        #we are going to print off a simple tab of the genes assigned to a KO for each organism
        OutHandle = open(self.OutputPath, 'w')
        #create header
        Header = "Accession\tKoala\tManual\n"
        OutHandle.write(Header)
        #now the matrix of annotations
        for ProteinObject in self.ProteinObjectsDictionary.values():
            PrintString = ProteinObject.GetLine()
            OutHandle.write("%s\n"%PrintString)
        OutHandle.close()
       
    def AddFiles(self, Koala, GeneToAccession, GeneToKO):
        #this is called by the iPython Notebook to mimic the ParseCommandLine() function
        #but from the notebook
        self.KoalaAnnotationPath = Koala
        self.ManualGeneAccessionPath = GeneToAccession
        self.ManualGeneKOPath = GeneToKO

    def SetIncludeAllAnnotations(self, Value):
        #this is called by the iPython Notebook to mimic the ParseCommandLine() function
        #but from the notebook
        #this specifically sets the variable that allows sub-par annotations to
        #be included. Default from the constructor is 0. meaning DO NOT include
        self.IncludeAllAssignments = Value


    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "k:m:n:a")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-m":
                self.ManualGeneAccessionPath = Value
            if Option == "-n":
                self.ManualGeneKOPath = Value 
            if Option == "-a":
                self.IncludeAllAssignments = 1
            if Option == "-k":
                self.KoalaAnnotationPath = Value
        # Error out, if we didn't see required options:
        
if __name__ == "__main__":
    Robot = CreatorClass()
    Robot.ParseCommandLine(sys.argv[1:])
    Robot.Main()
