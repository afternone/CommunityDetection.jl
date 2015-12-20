# type of undirected flow graph
type FlowGraph{V}
    graph::AbstractGraph{V}
    visit_prob::Vector{Float64} # nodes visit probability
    trans_prob::Vector{Float64} # edges transform probability
end

"construct undirected FlowGraph{V} from undirected AbstractGraph{V}"
function flow_graph{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)))
    FlowGraph(graph, visit_prob_undirected(graph, weights), trans_prob_undirected(graph, weights))
end

"calculate visit probability of nodes in undirected graphs"
function visit_prob_undirected{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)))
    @graph_requires graph vertex_list vertex_map adjacency_list edge_map

    !is_directed(graph) || error("graph must be undirected.")
    length(weights)==num_edges(graph) || error("length(weights) must equal num_edges(graph)")

    visit_prob = zeros(num_vertices(graph))
    for u in vertices(graph)
        for e in out_edges(u, graph)
            visit_prob[vertex_index(u, graph)] += weights[edge_index(e, graph)]
        end
    end
    visit_prob/2sum(weights)
end

"calculate transform probability of edges in undirected graphs"
function trans_prob_undirected{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)))
    !is_directed(graph) || error("graph must be undirected.")
    length(weights)==num_edges(graph) || error("length(weights) must equal num_edges(graph)")

    weights/2sum(weights)
end
