# using also http://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/


### * class GOparserError(Exception)
class GOparserError(Exception) :
    """Modified from the tutorial.pdf file of the Python documentation"""
    def __init__(self, value) :
        self.value = value
    def __str__(self) :
        return(repr(self.value))

### * class GOreader
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
                line = next(self.fileHandle).strip()
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
                        content = next(self.fileHandle).strip()
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

    def __next__(self) :
        return(self.next())

### * class GOnode
class GOnode :
    """GO tree node (GO term)"""

    def __init__(self, GOid, GOparents, GOchildren) :
        """Initialization.
        GOid is a string.
        GOparents and GOchildren are sets of strings.
        """
        self.GOid = GOid
        self.GOparents = GOparents
        self.GOchildren = GOchildren

### * class GOtree
class GOtree :
    """GO tree"""

    def __init__(self, listOfGOterms) :
        """Initialization
        listOfGOterms is something like:
        x = [x for x in GOreader("go.obo")]
        """
        self.GO = listOfGOterms
        self.checkNamespaces()
        self.checkUniquenessIds()
        self.checkNames()
        self.makeGOdict()
        self.buildTree()
        self.buildAltIdDict()
        
    def checkNamespaces(self) :
        """Check that each entry has only one namespace among the three allowed
        namespaces"""
        nsp = [x["namespace"] for x in self.GO]
        lnsp = [len(x) for x in nsp]
        lnsp = set(lnsp)
        if (not lnsp == set([1])) :
            raise GOparserError("GO term namespace error (len!=1)")
        nsp = [x[0] for x in nsp]
        nsp = set(nsp)
        if (not nsp == set(["cellular_component", "biological_process",
                            "molecular_function"])) :
            raise GOparserError(("Improper GO term namespaces: " +
                                 repr(list(nsp))))

    def checkUniquenessIds(self) :
        """Check that each entry has a unique id"""
        l = set([len(x["id"]) for x in self.GO])
        if (not l == set([1])) :
            raise GOparserError("Not all entries have exactly one id")
        ids = set([x["id"][0] for x in self.GO])
        if (not len(ids) == len(self.GO)) :
            raise GOparserError("Not all entries have a unique ids")

    def checkNames(self) :
        """Check that each entry has exactly one name"""
        l = set([len(x["name"]) for x in self.GO])
        if (not l == set([1])) :
            raise GOparserError("Not all entries have exactly one name")
        
    def makeGOdict(self) :
        """Make a dictionary (GOid: GOterm) from the GO list"""
        self.GOdict = dict(zip([x["id"][0] for x in self.GO], self.GO))
        if (not len(self.GOdict.keys()) == len(self.GO)) :
            raise GOparserError("Error while building the GO dictionary")

    def buildTree(self) :
        """Build the graph (tree) of GO entries"""
        self.GOtree = dict()
        for entry in self.GO :
            GOid = entry["id"][0]
            # get the previous record, if exists
            currentNode = self.GOtree.get(GOid, GOnode(GOid, set([]), set([])))
            # add the parents, if any
            for parent in entry.get("is_a", []) :
                parent = parent.split(" ! ", 1)[0]
                currentNode.GOparents.add(parent)
                # update the children sets of the parents
                parentNode = self.GOtree.get(parent,
                                             GOnode(parent, set([]), set([])))
                parentNode.GOchildren.add(GOid)
                self.GOtree[parent] = parentNode
            self.GOtree[GOid] = currentNode

    def buildAltIdDict(self) :
        """Build a dictionary (alt_id, real_id). All alt_ids are stored there,
        but for convenience the real_ids are also there as (real_id, real_id).
        This means that sending any id to the dictionary will return the 
        corresponding real id, even if the query itself was a real id."""
        altIdDict = dict()
        for GOid in self.GO :
            real_id = GOid["id"][0]
            alt_id = [real_id] + GOid.get("alt_id", [])
            for each in alt_id :
                assert altIdDict.get(each, "not_recorded") == "not_recorded"
                altIdDict[each] = real_id
        self.altIdDict = altIdDict

    def getAncestors(self, GOid, unique = False, depth = -1) :
        """Get the list of ancestors from a GO node (GO term)
        If unique is True, only returns unique ancestors.
        If depth = 0, doesn't return anything.
        If depth = n > 0, go as deep in the descendants (n recursive calls).
        If depth < 0, go all the way to the terminal leaves."""
        if (depth == 0) :
            return([])
        depth = depth - 1
        parents = list(self.GOtree.get(GOid, GOnode(set([]),
                                                    set([]),
                                                    set([]))).GOparents)
        ancestors = parents + []
        for p in parents :
            ancestors += self.getAncestors(p, unique, depth)
        if (unique) :
            return(list(set(ancestors)))
        else :
            return(ancestors)

    def getAncestorsGraph(self, GOid, depth = -1) :
        """Get the list of graph edges from a GO node to its ancestors
        If depth = 0, doesn't return anything.
        If depth = n > 0, go as deep in the descendants (n recursive calls).
        If depth < 0, go all the way to the terminal leaves."""
        if (depth == 0) :
            return([])
        depth = depth - 1
        parents = list(self.GOtree.get(GOid, GOnode(set([]),
                                                    set([]),
                                                    set([]))).GOparents)
        edges = [(GOid, x) for x in parents]
        for p in parents :
            edges += self.getAncestorsGraph(p, depth)
        return(list(set(edges)))

    def getGenealogyGraph(self, GOid, depthDescendants = -1, depthAncestors = -1) :
        """Get the graph of both descendants and ancestors of a GO node.
        Depth can be set for descendants and ancestors.
        A negative depth means to go all the way to the terminal nodes."""
        g = self.getDescendantsGraph(GOid, depthDescendants)
        g += self.getAncestorsGraph(GOid, depthAncestors)
        return(g)

    def showGOnames(self, listOfGOid) :
        """Get the names of GO terms"""
        dicts = [self.GOdict.get(x, dict(name = ["NA"])) for x in listOfGOid]
        return([x["name"][0] for x in dicts])

    def showGOnamesGraph(self, graphOfGOid) :
        """Get the names of GO terms in the output from 
        self.getAncestorsGraph"""
        o = [(self.GOdict[x[0]]["name"][0],
              self.GOdict[x[1]]["name"][0]) for x in graphOfGOid]
        return(o)

    def graphvizGraph(self, graphOfGOid, names = False) :
        """Produce a graphiz input file from the output of 
        self.getAncestorsGraph
        If names if True, add the corresponding names"""
        if (names) :
            graphData = list()
            graphOfGOidNames = self.showGOnamesGraph(graphOfGOid)
            for (l1, l2) in zip(graphOfGOid, graphOfGOidNames) :
                graphData.append((l1[0] + "\\n" + l2[0],
                                  l1[1] + "\\n" + l2[1]))
            graphOfGOid = graphData
        labels = dict()
        n_labels = 0
        for l in graphOfGOid :
            for ll in l :
                if (ll in labels.keys()) :
                    pass
                else :
                    labels[ll] = "label" + str(n_labels)
                    n_labels += 1
        o = "digraph G {\n"
        o += "    rankdir=BT\n"
        for l in graphOfGOid :
            o += "    " + labels[l[0]] + " -> " + labels[l[1]] + ";\n"
        for (k, v) in labels.items() :
            o += "    " + v + " [shape=box, label=\"" + k + "\"];\n" 
        o += "}\n"
        return(o)            

    def getDescendants(self, GOid, unique = False, depth = -1) :
        """Get the list of descendants from a GO node (GO term)
        If unique is True, only returns unique decendants.
        If depth = 0, doesn't return anything.
        If depth = n > 0, go as deep in the descendants (n recursive calls).
        If depth < 0, go all the way to the terminal leaves."""
        if (depth == 0) :
            return([])
        depth = depth - 1
        children = list(self.GOtree.get(GOid, GOnode(set([]),
                                                    set([]),
                                                     set([]))).GOchildren)
        descendants = children + []
        for c in children :
            descendants += self.getDescendants(c, unique, depth)
        if (unique) :
            return(list(set(descendants)))
        else :
            return(descendants)

    def getDescendantsGraph(self, GOid, depth = -1) :
        """Get the list of graph edges from a GO node to its descendants
        If depth = 0, doesn't return anything.
        If depth = n > 0, go as deep in the descendants (n recursive calls).
        If depth < 0, go all the way to the terminal leaves."""
        if (depth == 0) :
            return([])
        children = list(self.GOtree.get(GOid, GOnode(set([]),
                                                    set([]),
                                                        set([]))).GOchildren)
        depth = depth - 1
        edges = [(x, GOid) for x in children]
        for c in children :
            edges += self.getDescendantsGraph(c, depth)
        return(list(set(edges)))

    def showGraph(self, graph) :
        """Display a graph. Relies on dot and the Image module on Linux.
        Not tested on Windows."""
        import uuid
        temp = str(uuid.uuid4())
        f = open(temp + ".dot", "w")
        f.write(self.graphvizGraph(graph, True))
        f.close()
        import subprocess
        c = ["dot", "-Tpng", temp + ".dot", "-o", temp + ".png"]
        p = subprocess.Popen(c)
        p.wait()
        import Image
        Image.open(temp + ".png").show()
        import os
        os.remove(temp + ".dot")
        os.remove(temp + ".png")

### * class GOroot
class GOroot :
    """A simple class to hold the roots of the tree"""

    def __init__(self) :
        self.BP = ["GO:0008150", "GO:0007582", "GO:0000004"]
        self.CC = ["GO:0005575", "GO:0008372"]
        self.MF = ["GO:0003674", "GO:0005554"]

# test

import cProfile

def count(i) :
    o = dict()
    for x in i :
        o[x] = o.get(x, 0)
        o[x] += 1
    return(o)

#a = GOreader("go.obo")

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

#GO = [x for x in a]

a = GOtree([x for x in GOreader("go.obo")])
goid = "GO:0046405"
g = a.graphvizGraph(a.getAncestorsGraph(goid), True)
f = open("toto.dot", "w")
f.write(g)
f.close()

# dot toto.dot -Tpdf > toto.pdf

# test for extracting parents
#import random
#n = 100000
#l = list(a.GOdict.keys())
#GOs = [random.choice(l) for x in xrange(n)]
#cProfile.run('x = [a.getAncestors(x) for x in GOs]')
#          21742265 function calls (14594844 primitive calls) in 18.963 seconds

#    Ordered by: standard name

#    ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#         1    0.069    0.069   18.963   18.963 <string>:1(<module>)
# 7247421/100000   16.322    0.000   18.895    0.000 GO_parser.py:148(getAncestors)
#   7247421    1.813    0.000    1.813    0.000 GO_parser.py:70(__init__)
#         1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
#   7247421    0.759    0.000    0.759    0.000 {method 'get' of 'dict' objects}
