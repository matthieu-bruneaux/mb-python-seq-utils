<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#sec-1">1. Gene Ontology tools</a>
<ul>
<li><a href="#sec-1-1">1.1. Classes defined in <b>GO<sub>parser</sub>.py</b></a></li>
<li><a href="#sec-1-2">1.2. <b>GOtree</b> attributes of interest</a></li>
<li><a href="#sec-1-3">1.3. <b>GOtree</b> methods of interest</a></li>
<li><a href="#sec-1-4">1.4. TO DO</a></li>
</ul>
</li>
<li><a href="#sec-2">2. Fasta sequence file tools</a></li>
</ul>
</div>
</div>
A set of various short scripts written in Python and used for sequence 
manipulation, gene ontology manipulation and such.

If you use those scripts and find any error, please report it!

# Gene Ontology tools

Those tools are located in the **gene<sub>ontology</sub>** folder. The most interesting
file is **GO<sub>parser</sub>.py**: for now there is a class implementing a simple GO obo
file format reader and a GOtree class with methods to retrieve all ancestors
from a given GO id.

GO tree is built based on the *is<sub>a</sub>* fields, and no checking against cyclic
graph is performed. The resulting graph ([example](https://github.com/matthieu-bruneaux/python-bioinformatic-utils/raw/master/gene_ontology/toto.pdf)) can be produced using
**graphviz**

In terms of performances, a GO obo file with 41865 entries (**go.obo**) can be
loaded in a list using **GOreader** in 1.46s on my computer. After building the
tree, all ancestors can be retrieved for a random set of 100000 (1e5) GO ids in
19s (linear run time).

A repository which looks very interesting and with much more advanced code that
this repository is [goatools](https://github.com/tanghaibao/goatools).

## Classes defined in **GO<sub>parser</sub>.py**

-   **`GOparserError`:** error raised when a problem occurred in the parsing of the
    GO data

-   **`GOreader`:** parser for files in GO obo format. It is an iterable, so list
    comprehensions can be used to load all the data from a GO obo
    file into a list.

-   **`GOnode`:** structure used in GOtree to represent one node in the GO
    hierarchy

-   **`GOtree`:** main class in the module. When initialized, an instance takes a
    list of GO terms as can be produced with GOreader and a list
    comprehension.

-   **`GOroot`:** a simple class to hold the root ids of the tree

## **GOtree** attributes of interest

-   **`GO`:** the original GO term list used to build the tree

-   **`GOdict`:** a dictionary mapping each GO `id` to the full description of the
    term

-   **`altIdDict`:** a dictionary mapping each GO `id` and `alt_id` to the main
    `id` they refer to. When analyzing a list of GO ids which
    might contain ids which are stored as `alt_id` in the obo
    file, this dictionary can be used to convert any id value to
    the currently accepted `id` for the GO term

## **GOtree** methods of interest

-   **`getAncestors(GOid, unique = FALSE, depth = -1)`:** return a list of the GO
    ids of the ancestors of `GOid`

-   **`getDescendants(GOid, unique = FALSE, depth = -1)`:** return a list of the
    descendants of `GOid`

-   **`showGOnames(listOfGOid)`:** return a list of GO term names corresponding to
    the GO ids given as input

## TO DO

### TODO Add some example scripts

### TODO Generate a proper documentation

# Fasta sequence file tools

Those tools are in the **sequences** folder. Each script performs a very simple
task.
