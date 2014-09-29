Python bioinformatic utilities
==============================

A set of various short scripts written in Python and used for sequence 
manipulation, gene ontology manipulation and such.

If you use those scripts and find any error, please report it!

Gene Ontology tools
-------------------

Those tools are located in the **gene_ontology** folder. The most interesting
file is **GO_parser.py**: for now there is a class implementing a simple GO obo
file format reader and a GOtree class with methods to retrieve all ancestors
from a given GO id.

GO tree is built based on the *is_a* fields, and no checking against cyclic
graph is performed. The resulting graph
([example](https://github.com/matthieu-bruneaux/python-bioinformatic-utils/raw/master/gene_ontology/toto.pdf))
can be produced using **graphviz**.

A repository which looks very interesting and with much more advanced code that
this repository is [goatools](https://github.com/tanghaibao/goatools).

Fasta sequence file tools
-------------------------

Those tools are in the **sequences** folder. Each script performs a very simple
task.
