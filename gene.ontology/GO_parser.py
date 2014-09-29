# using also http://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/

class GOreader :
    """GO obo file reader"""

    def __init__(self, filename) :
        """Initialization"""
        self.filename = filename
        self.fileHandle = open(filename, mode = "r")
        self.termFound = False
        self.inTerm = False

    def __iter__(self) :
        return(self)

    def next(self) :
        """Get the next GO term"""
        while (not self.termFound) :
            line = self.fileHandle.next().strip()
            if (line == "[Term]") :
                self.termFound = True
        self.inTerm = True
        term = dict()
        line = self.fileHandle.next().strip()
        if (line == "[Term]" or line == "[Typedef]") :
            self.inTerm = False
            if (line == "[Term]") :
                self.inTerm = True
        while (self.inTerm) :
            if (line == "[Term]" or line == "[Typedef]") :
                self.inTerm = False
            elif (line != "") :    
                (a, b) = line.split(": ", 1)
                # add error if split doesn't work (length != 2)
                term[a] = term.get(a, [])
                term[a].append(b)
                line = self.fileHandle.next().strip()
            else :
                line = self.fileHandle.next().strip()
        return(term)
        
a = GOreader("goslim_generic.obo")


