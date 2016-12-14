# ------------------------------------------------------------------------------------
# ------------------------- Legend for the EDP benchmarks    -------------------------
# ------------------------------------------------------------------------------------

The BP algorithm reference paper can be found at:  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0145222
For a beta version of the code please contact: cdebacco@santafe.edu
Eventually it will be uploaded to github...

This dataset can be found at:
http://becool.info.ucl.ac.be/lsgraph

under "lsgraph_bench_edp.zip". Thanks to Pham Quang Dung.

The titles denote the graph type, and "ins1" or "ins2" denote the instance number: two different instances have the same graph but different set of (S,R) pairs.

Row 1: $1=N=# of vertices; $2=E=# of edges

All the rows with a third column are the edges, the third column denotes a weight. But in all these examples is set by default to 0, because they all consider unweighted networks.

Example:
202 217 0 
Denotes an edge from node 202 to node 217.

The last rows of the file, only 2 columns, are the (Source, Receiver) pairs of the instance. These should be fixed and one should look for the edge disjoint paths connecting these fixed pairs.
Example:
78 9
denotes source node=78 and receiver node=9

If you want to compare with the results reported in [Altarelli, Fabrizio, Alfredo Braunstein, Luca Dall’Asta, Caterina De Bacco, and Silvio Franz. "The edge-disjoint path problem on random graphs by message-passing." PloS one 10, no. 12 (2015): e0145222] you should take the max,min and average number of edge-disjoint path accomodated (among the fixed (S,R) pairs) over 20 runs of your algorithm for each instance ("ins1" or "ins2").

If you use this benchmark please cite:

1. [Altarelli, Fabrizio, Alfredo Braunstein, Luca Dall’Asta, Caterina De Bacco, and Silvio Franz. "The edge-disjoint path problem on random graphs by message-passing." PloS one 10, no. 12 (2015): e0145222]
url: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0145222

2.[Pham, Quang Dung, Yves Deville, and Pascal Van Hentenryck. "LS (Graph): a constraint-based local search for constraint optimization on trees and paths." Constraints 17, no. 4 (2012): 357-408]
