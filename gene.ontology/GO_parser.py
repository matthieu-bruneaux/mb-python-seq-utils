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

# test

import cProfile

a = GOreader("go.obo")

# cProfile.run('x = [x for x in a]')
#       1989406 function calls in 1.461 seconds

# Ordered by: standard name

# ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#      1    0.021    0.021    1.461    1.461 <string>:1(<module>)
#      1    0.000    0.000    0.000    0.000 GO_parser.py:12(__iter__)
#  41866    1.110    0.000    1.439    0.000 GO_parser.py:15(next)
# 465926    0.050    0.000    0.050    0.000 {method 'append' of 'list' objects}
#      1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
# 465926    0.056    0.000    0.056    0.000 {method 'get' of 'dict' objects}
# 465926    0.151    0.000    0.151    0.000 {method 'split' of 'str' objects}
# 549759    0.072    0.000    0.072    0.000 {method 'strip' of 'str' objects}

GO = [x for x in a]
