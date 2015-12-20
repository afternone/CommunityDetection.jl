# CommunityDetection

[![Build Status](https://travis-ci.org/afternone/CommunityDetection.jl.svg?branch=master)](https://travis-ci.org/afternone/CommunityDetection.jl)

INTRODUCTION
============
Community detection in Julia.
This package inspired by [louvain-igraph](https://github.com/vtraag/louvain-igraph).
It relies on `Graphs.jl` for it to function. Besides the relative
flexibility of the implementation, it also scales well, and can be run on graphs
of millions of nodes (as long as they can fit in memory). 


The core function is 
``optimize_partition`` which finds the optimal partition using the louvain algorithm
for a number of different methods. The methods currently implemented are:

* Modularity.
  This method compares the actual graph to the expected graph, taking into
  account the degree of the nodes [1]. The expected graph is based on a
  configuration null-model. Notice that we use the non-normalized version (i.e.
  we don't divide by the number of edges), so that this Modularity values
  generally does not fall between 0 and 1. The formal definition is

  ```
  H = sum_ij (A_ij - k_i k_j / 2m) d(s_i, s_j),
  ```

  where `A_ij = 1` if there is an edge between node `i` and `j`, `k_i` is the degree of
  node `i` and `s_i` is the community of node i.

* RBConfiguration.
  This is an extension of modularity which includes a resolution parameter [2].
  In general, a higher resolution parameter will lead to smaller communities.
  The formal definition is

  ```
  H = sum_ij (A_ij - gamma k_i k_j / 2m) d(s_i, s_j),
  ```

  where `gamma` is the resolution value, and the other variables are the same as
  for Modularity.

* RBER.
  A variant of the previous method that instead of a configuration null-model
  uses a Erdös-Rényi null-model in which each edge has the same probability of
  appearing [2]. The formal definition is

  ```
  H = sum_ij (A_ij - gamma p) d(s_i, s_j),
  ```

  where `p` is the density of the graph, and the other variables are the same as
  for Modularity, with `gamma` a resolution parameter.


* CPM.
  This method compares to a fixed resolution parameter, so that it finds
  communities that have an internal density higher than the resolution
  parameter, and is separated from other communities with a density lower than
  the resolution parameter [3].The formal definition is

  ```
  H = sum_ij (A_ij - gamma ) d(s_i, s_j),
  ```

  with `gamma` a resolution parameter, and the other variables are the same as for
  Modularity.

* Significance.
  This is a probabilistic method based on the idea of assessing the probability
  of finding such dense subgraphs in an (ER) random graph [4]. The formal
  definition is

  ```
  H = sum_c M_c D(p_c || p)
  ```

  where `M_c` is the number of possible edges in community `c`, i.e. `n_c (n_c - 1)/2`
  for undirected graphs and twice that for directed grahs with `n_c` the size of
  community `c`, `p_c` is the density of the community `c`, and `p` the general density
  of the graph, and `D(x || y)` is the binary Kullback-Leibler divergence.

* Surprise.
  Another probabilistic method, but rather than the probability of finding dense
  subgraphs, it focuses on the probability of so many edges within communities
  [5, 6]. The formal definition is

  ```
  H = m D(q || <q>)
  ```

  where `m` is the number of edges, `q` is the proportion of edges within
  communities (i.e. `sum_c m_c / m`) and `<q>` is the expected proportion of edges
  within communities in an Erdős–Rényi graph.

* Infomap (Map Equation)
A method based on infomation theory, 
you can explore the mechanics of the map equation [in here](http://www.mapequation.org/)

This package also implement label propagation algorithm and neighbor strength driven label propagation algorithm.

INSTALLATION
============

```
julia> Pkg.clone("git://github.com/afternone/CommunityDetection.jl.git")
```

USAGE
=====

To start, make sure to import the packages:
```
using CommunityDetection
using GraphPlot # for testing graph
```

We'll create a simple graph for testing purposes:
```
g = graphfamous("karate")
```

For simply finding a partition use:
```
# construct modularity partition
mp = mpartition(g)
optimize_partition(mp)
```

In case you want to use a weighted graph
```
edge_weights = ones(num_edge(g))
mp = mpartition(g, edge_weights)
optimize_partition(mp)
```
Please note that not all methods are necessarily capable of handling weighted
graphs.

Notice that ``mp`` now contains the membership of nodes
```
mp.membership
```

You can also find partition using infomap algorithm
```
# construct flow partition
fp = flow_partition(g)
optimize_partition(fp)
```

Some other examples
```
# Significance method
mp = mpartition(g)
optimize_partition(mp, method = :Significance)

# label propagation
membership1 = lpa(g)

# neighbor strength driven label propagation
membership2 = nsdlpa(g)

# calculate modularity
modularity(g, membership1)

# calculate NMI of two partition
nmi(membership1, membership2)

# calculate VOI of two partition
voi(membership1, membership2)
```


REFERENCES
==========

Please cite the references appropriately in case they are used.

1. Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. Fast unfolding
   of communities in large networks. J. Stat. Mech. 2008, P10008 (2008).
2. Newman, M. & Girvan, M. Finding and evaluating community structure in networks.
   Physical Review E 69, 026113 (2004).
3. Reichardt, J. & Bornholdt, S. Partitioning and modularity of graphs with arbitrary
   degree distribution. Physical Review E 76, 015102 (2007).
4. Traag, V. A., Van Dooren, P. & Nesterov, Y. Narrow scope for resolution-limit-free
   community detection. Physical Review E 84, 016114 (2011).
5. Traag, V. A., Krings, G. & Van Dooren, P. Significant scales in community structure.
   Scientific Reports 3, 2930 (2013).
6. Aldecoa, R. & Marín, I. Surprise maximization reveals the community structure
   of complex networks. Scientific reports 3, 1060 (2013).
7. Traag, V.A., Aldecoa, R. & Delvenne, J.-C. Detecting communities using Asymptotical
   Surprise. Forthcoming (2015).
8. Mucha, P. J., Richardson, T., Macon, K., Porter, M. A. & Onnela, J.-P.
   Community structure in time-dependent, multiscale, and multiplex networks.
   Science 328, 876–8 (2010).

