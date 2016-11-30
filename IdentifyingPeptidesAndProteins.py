#parse PSM results to get peptide and protein identifications


#this is written in Python 3 

HelpOutput = """
python IdentifyingPeptidesAndProteins.py -d /path/to/files -o file.txt

Required Parameters:
-d   <dir>  This is a path to the directory that holds all peptide identifications
-o   <file> This is a file that tells which PSM files belong to a specific organism

Example: python .\IdentifyingPeptidesAndProteins.py -d ..\MSGF.Searches -o FileToOrganisms.txt
"""


import os
import sys
import getopt
import string



class Protein:
    def __init__(self, accession):
        self.ProteinAccession = accession
        self.KeggOrthologAccession = ""
        self.Peptides = {} #key = sequence, value = spectrum count

    def AddPeptide(self, Peptide):
        if not Peptide in self.Peptides:
            self.Peptides[Peptide] = 0
        self.Peptides[Peptide] += 1


class Organism:
    def __init__(self, Name):
        self.OrganismName = Name
        self.Proteins = {} #key = Accession, value= object

    def AddPeptide(self, Peptide, ProteinAcc):
        if not ProteinAcc in self.Proteins:
            self.Proteins[ProteinAcc] = Protein(ProteinAcc)
        self.Proteins[ProteinAcc].AddPeptide(Peptide)
    def GetProteinCount(self):
        #simply return the number of proteins
        return len(self.Proteins)

class ParserClass:

    def __init__(self):
        "nothing to put in really"
        self.DirectoryOfPSMs = "" # this dir holds all the PSM results (txt files)
        self.FileToOrganismsPath  = ""
        self.QvalueCutoff = 0.001
        self.OrganismObjectsDictionary = {} #key = Name, value= object
        self.FileToOrganismDictionary = {} #key = PSM file stub, value = organism
        self.OutputPath = "renameYourOutput.txt"

    def GetOrganismObjects(self):
        return self.OrganismObjectsDictionary

    def SetQvalue(self, qvalue):
        #for the ipython notebook
        self.QvalueCutoff = qvalue

    def Main(self):
        #0. Figure out which PSM files belong with which organisms
        self.ParseFileToOrganism(self.FileToOrganismsPath)
        #1. parse all the PSM files just putting peptides into proteins
        ItemsInDir = os.listdir(self.DirectoryOfPSMs)
        for Item in ItemsInDir:
            #I am expecting .txt files
            if not Item[-4:] == ".txt":
                continue
            Path = os.path.join(self.DirectoryOfPSMs, Item)
            FileStub = os.path.splitext(Item)[0] #gets the root filename, and removes the extension
            #have to do some more cleanup to get back a good file name
            FileStub = FileStub.replace("_msgfdb_fht", "")
            if not FileStub in self.FileToOrganismDictionary:
                print ("Can't associate this file %s with an organism."%Item)
                print ("In directory %s, but not listed in association file %s"%(self.DirectoryOfPSMs, self.FileToOrganismsPath))
                continue
            OrganismName = self.FileToOrganismDictionary[FileStub]
            if not os.path.isfile(Path):
                continue
            self.ParsePSMFile(Path, OrganismName)

    def ParseFileToOrganism(self, Path):
        #simple populating the FileToOrganismDictionary variable so that we can use this knowledge later on
        Handle = open(Path, 'r')
        print ("I got this path %s"%Path)
        Header = Handle.readline() #pop it off because I don't need it
        FileCounter = 0
        OrganismCounter = 0
        for Line in Handle:
            #[0] = organism name
            #[1] = comma separated list of file stubs
            Bits = Line.strip().split("\t")
            #make the organism object
            OrganismName = Bits[0]
            if not OrganismName in self.OrganismObjectsDictionary:
                self.OrganismObjectsDictionary[OrganismName] = Organism(OrganismName)
                OrganismCounter += 1
            #parse out the files
            Files = Bits[1].split(",")
            for File in Files:
                #probably have to strip out some white space here
                File = File.strip()
                self.FileToOrganismDictionary[File] = OrganismName
                FileCounter += 1
        print("Associating %s files from %s organisms"%(FileCounter, OrganismCounter))

    def ParsePSMFile(self, Path, OrganismName):

        Org = self.OrganismObjectsDictionary[OrganismName]
        Handle = open(Path, 'r')
        Header = Handle.readline()
        for Line in Handle:
            Bits = Line.strip().split("\t")
            # [4] is charge
            #[9] is the peptide string
            #[10] is protein
            #[14] is SPecEvalue
            #[17] is q value
            peptide = Bits[9][2:-2] #this last bit is to strip out the prefix/suffix
            qvalue = float(Bits[17])
            protein = Bits[10]
            #first a hard filter on the qvalue. We simply don't bother to look at crap
            if qvalue > self.QvalueCutoff:
                continue
            #now we don't want to include all the contaminants or decoys
            if protein[:3] in ["Con", "XXX"]:
                continue
            ### Currently adding the entire protein fasta accession set
            ## sp|P54547|G6PD_BACSU
            ## ref|YP_001233830.1
            ## tr|J0JN52|J0JN52_ALCFA
            Org.AddPeptide(peptide, protein)


       
    def AddFiles(self, PSM_dir, FileToOrganisms):
        #this is called by the iPython Notebook to mimic the ParseCommandLine() function
        #but from the notebook
        self.FileToOrganismsPath = FileToOrganisms
        self.DirectoryOfPSMs = PSM_dir


    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:o:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                self.DirectoryOfPSMs = Value
            if Option == "-o":
                self.FileToOrganismsPath = Value 
        
if __name__ == "__main__":
    Robot = ParserClass()
    Robot.ParseCommandLine(sys.argv[1:])
    Robot.Main()
