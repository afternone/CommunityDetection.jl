# type of directed flow graph
type DiFlowGraph{V}
    graph::AbstractGraph{V}
    tau::Float64 # damping parameter
    visit_prob::Vector{Float64} # nodes visit probability
    trans_prob::Vector{Float64} # edges transform probability
end

"construct directed DiFlowGraph{V} from directed AbstractGraph{V}"
function diflow_graph{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)); τ=0.15)
    trans_prob = trans_prob_directed(graph, weights)
    visit_prob = visit_prob_directed(graph, trans_prob, τ=τ, ϵ=sqrt(eps()), N=1000)
    DiFlowGraph(graph, τ, visit_prob, trans_prob)
end

function _diflow_graph{V,T<:Real}(graph::AbstractGraph{V}, trans_prob::Vector{T}; τ=0.15)
    visit_prob = visit_prob_directed(graph, trans_prob, τ=τ, ϵ=sqrt(eps()), N=1000)
    DiFlowGraph(graph, τ, visit_prob, trans_prob)
end

"calculate relative transform probability of edges in a directed graphs"
function trans_prob_directed{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)))
    @graph_requires graph vertex_list vertex_map adjacency_list edge_map

    is_directed(graph) || error("graph must be directed.")
    length(weights)==num_edges(graph) || error("length(weights) must equal num_edges(graph)")

    trans_prob = zeros(num_edges(graph))

    for u in vertices(graph)
        sumw = 0.0
        for e in out_edges(u, graph)
            sumw += weights[edge_index(e, graph)]
        end
        if sumw > 0.0
            for e in out_edges(u, graph)
                e_idx = edge_index(e, graph)
                trans_prob[e_idx] = weights[e_idx]/sumw
            end
        end
    end

    trans_prob
end

"""
calculate visit probability of each node in directed networks
`tau` is damping parameter,  `epsilon` is convergence precision, `N` is maximum iterations
"""
function visit_prob_directed{V,T<:Real}(graph::AbstractGraph{V}, trans_prob::Vector{T}; τ=0.15, ϵ=10e-5, N=1000)
    @graph_requires graph vertex_list vertex_map adjacency_list edge_map

    is_directed(graph) || error("graph must be directed.")
    length(trans_prob)==num_edges(graph) || error("length(relative_trans_prob) must equal num_edges(graph)")

    n = num_vertices(graph)
    p = fill(1/n, n) # initial visit probability
    p1 = zeros(n)
    i = 0 # initial iterations
    δp = 1.0 # δp > ϵ for into while loop at first time

    while δp > ϵ && i < N
        i += 1
        dp = 0.0
        for u in vertices(graph)
            if out_degree(u, graph) == 0
                dp += (1.0 - τ) * p[vertex_index(u, graph)] / n
            end
        end
        for u in vertices(graph)
            u_idx = vertex_index(u, graph)
            p1[u_idx] = dp + τ/n
            for e in in_edges(u, graph)
                v = source(e, graph)
                v_idx = vertex_index(v, graph)
                e_idx = edge_index(e, graph)
                p1[u_idx] += (1.0-τ) * trans_prob[e_idx] * p[v_idx]
            end
        end
        δp = sumabs(p1-p)
        p[:] = p1[:]
    end

    p
end
