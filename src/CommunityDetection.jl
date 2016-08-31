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
	nmi, nvoi, nminvoi, rnmi, modularity, conductance, min_conductance


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

end # module
