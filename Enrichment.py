import scipy.stats as stats


class PathwayEnrichment:
    #mostly a convenience class for making plots. Could all be done with a set
    #of embedded dictionaries. But no one likes those.
    def __init__(self, PathwayID, PathwayDescription):
        self.PathwayID = PathwayID # the kegg identifier
        self.PathwayDescription = PathwayDescription #more of 'glycolysis'
        self.EnrichmentByOrg = {} #key = organism name, value = pvalue
    def GetEnrichmentList(self):
        #return just a simple list of the pvalues, stripped from the association with an organism
        return list(self.EnrichmentByOrg.values())
    def AddValue(self, org, pvalue):
        self.EnrichmentByOrg[org] = pvalue
    def GetMedian(self):
        #this function looks at the pvalues and finds the median one.
        ListOfPvalues = self.GetEnrichmentList()
        ListOfPvalues.sort() #sort in place
        Length = len(ListOfPvalues)
        HalfWay = int(Length / 2 ) #probably not truly accurate because I don't average the middle if even length, 
        return ListOfPvalues[HalfWay]


