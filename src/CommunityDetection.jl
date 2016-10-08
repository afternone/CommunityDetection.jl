module CommunityDetection

using Graphs
import Graphs.graph
using StatsBase

export
    AbstractPartition,
    flow_partition,
    diflow_partition,
    mpartition,
    optimize_partition!,
    find_partition!,
    multi_greedy!,
    lpa, nsdlpa, hlpa, hlpa_record, slpa, getgrp,
	nmi, nvoi, nminvoi, rnmi, modularity, conductance, min_conductance,
  internal_density, edges_inside, average_degree, fomd,
  triangle_participation, triangle_participation_ratio, triangle_participation_ratio,
  expansion, cut_ratio, normalized_cut, odf, max_odf, average_odf,
  flake_odf, separability, density,
  prob_metric_cluster, prob_metric_graph, pair_count, rand_index,
  mirkin_metric, jaccard_index, adjusted_rand_index,
  f1_score, fbeta_score, recall, classification_report,
  confusion_matrix, jaccard_similarity, hamming_loss, accuracy,
  adjusted_mutual_info_score, mutual_info_score, normalized_mutual_info_score,
  modularity_density


include("AbstractPartition.jl")
include("FlowGraph.jl")
include("FlowPartition.jl")
include("DiFlowGraph.jl")
include("DiFlowPartition.jl")
include("MGraph.jl")
include("MPartition.jl")
include("utils.jl")
include("nmi.jl")
include("metrics.jl")
include("Optimiser.jl")

include("modularity.jl")

include("label_propagation.jl")
include("mlpa.jl")
include("slpa.jl")
include("MultiGreedy.jl")
include("modularity_density.jl")

end # module
