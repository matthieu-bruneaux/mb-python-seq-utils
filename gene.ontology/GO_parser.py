# using also http://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/

class GOreader :
    """GO obo file reader"""

    def __init__(self, filename) :
        """Initialization"""
        self.filename = filename
        self.fileHandle = open(filename, mode = "r")
        self.newTerm = False

    def __iter__(self) :
        return(self)

    def next(self) :
        """Get the next GO term"""
        while (True) :
            if (not self.newTerm) :
                # read the next line if we are not already at the beginning
                # of a new term
                line = self.fileHandle.next().strip()
            else :
                # we know from previous calls that we are at the beginning
                # of a new term
                line = "[Term]"
            if (line == "[Term]") :
                # new term
                self.newTerm = False
                term = dict()
                while (True) :
                    try :
                        content = self.fileHandle.next().strip()
                        if (content != "[Term]" and content != "[Typedef]") :
                            if (content != "") :
                                # new entry
                                (k, v) = content.split(": ", 1)
                                term[k] = term.get(k, [])
                                term[k].append(v)
                            else :
                                # empty entry
                                pass
                        elif (content == "[Term]") :
                            # new term
                            self.newTerm = True
                            # end of term
                            return(term)
                        else :
                            # simple end of term, next line is not a new term
                            return(term)
                    except StopIteration :
                        # end of file while in a term
                        return(term)
            else :
                # not a new term
                pass
        
a = GOreader("goslim_generic.obo")


