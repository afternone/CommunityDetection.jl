type MGraph{V}
    graph::AbstractGraph{V}
    node_sizes::Vector{Int}
    edge_weights::Vector{Float64}
    strength_in::Vector{Float64}
    strength_out::Vector{Float64}
    node_self_weights::Vector{Float64}
    total_weight::Float64
    total_size::Int
    correct_self_loops::Bool
    density::Float64
end


function mgraph{V,T}(g::AbstractGraph{V}, edge_weights::Vector{T}=ones(num_edges(g)),
                     node_sizes::Vector{Int}=ones(Int,num_vertices(g)), correct_self_loops::Bool=false)
    length(edge_weights) == num_edges(g) || error("wrong edge_weights length")
    length(node_sizes) == num_vertices(g) || error("wrong node_sizes length")

    m = num_edges(g)

    # determine total weight in the graph
    total_weight = 0.0
    for e in edges(g)
        e_idx = edge_index(e,g)
        total_weight += edge_weights[e_idx]
    end

    n = num_vertices(g)


    total_size = 0
    strength_in = Array(T, n)
    strength_out = Array(T, n)

    for u in vertices(g)
        u_idx = vertex_index(u,g)

        # determine total size in the graph
        total_size += node_sizes[u_idx]

        # determine strength in
        strength = 0.0
        for e in in_edges(u,g)
            e_idx = edge_index(e,g)
            strength += edge_weights[e_idx]
        end
        strength_in[u_idx] = strength

        # determine strength out
        strength = 0.0
        for e in out_edges(u,g)
            e_idx = edge_index(e,g)
            strength += edge_weights[e_idx]
        end
        strength_out[u_idx] = strength
    end

    # calculate density
    normalise = correct_self_loops ? total_size*total_size : total_size*(total_size-1)
    density = is_directed(g) ? total_weight/normalise : 2*total_weight/normalise

    # Set default self_weights of the total weight of any possible self-loops
    node_self_weights = zeros(T, n)
    for e in edges(g)
        e_idx = edge_index(e,g)
        u = source(e,g)
        v = target(e,g)
        u_idx = vertex_index(u,g)
        v_idx = vertex_index(v,g)
        if u_idx == v_idx
            node_self_weights[u_idx] = edge_weights[e_idx]
        end
    end

    MGraph(g, node_sizes, edge_weights, strength_in, strength_out, node_self_weights,
           total_weight, total_size, correct_self_loops, density)
end

# weight of vertex to community
function weight_to_community{V}(mg::MGraph{V}, u::V, comm::Int, membership::Vector{Int})
    g = mg.graph
    total_w = 0.0

    for e in out_edges(u,g)
        e_idx = edge_index(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        if membership[v_idx] == comm
            # Self loops appear twice here if the graph is undirected, so divide by 2.0 in that case
            w = mg.edge_weights[e_idx]
            if u == v && !is_directed(g)
                w /= 2
            end
            total_w += w
        end
    end

    total_w
end

# weight of vertex from community
function weight_from_community{V}(mg::MGraph{V}, u::V, comm::Int, membership::Vector{Int})
    g = mg.graph
    total_w = 0.0

    for e in in_edges(u,g)
        e_idx = edge_index(e,g)
        v = source(e,g)
        v_idx = vertex_index(v,g)
        if membership[v_idx] == comm
            # Self loops appear twice here if the graph is undirected, so divide by 2.0 in that case
            w = mg.edge_weights[e_idx]
            if u == v && !is_directed(g)
                w /= 2
            end
            total_w += w
        end
    end

    total_w
end
