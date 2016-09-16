function hlpa{V}(g::AbstractGraph{V}; lambda::Float64=1.0)
    collapsed_edge_weights = Dict{Int,Float64}[]
    sizehint!(collapsed_edge_weights, num_vertices(g))
    collapsed_weights = Float64[]
    sizehint!(collapsed_weights, num_edges(g))
    s = similarity(g)
    m = label_propagation(g, s, lambda)
    nb_comm1 = maximum(m)
    h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
    m1 = label_propagation(h, collapsed_weights, lambda)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = label_propagation(h, collapsed_weights, lambda)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = label_propagation(h, collapsed_weights, lambda)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = label_propagation(h, collapsed_weights, lambda)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = label_propagation(h, collapsed_weights, lambda)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m
end

"""
Community detection using the label propagation algorithm (see [Raghavan et al.](http://arxiv.org/abs/0709.2938)).
`g`: input Graph
`maxiter`: maximum number of iterations
return : vertex assignments and the convergence history
"""
function label_propagation(g, weights, lambda; maxiter=1000)
    n = num_vertices(g)
    label = collect(1:n)
    active_nodes = IntSet(vertices(g))
    c = NeighComm(collect(1:n), fill(-1.,n), 1)
    random_order = Array(Int, n)
    i = 0
    while !isempty(active_nodes) && i < maxiter
      num_active = length(active_nodes)
        i += 1

        # processing nodes in random order
        for (j,node) in enumerate(active_nodes)
            random_order[j] = node
        end
        range_shuffle!(1:num_active, random_order)
        @inbounds for j=1:num_active
            u = random_order[j]
            old_comm = label[u]
            label[u] = vote!(g, weights, label, c, u, lambda)
            if old_comm != label[u]
                for v in out_neighbors(u, g)
                    push!(active_nodes, v)
                end
            else
                delete!(active_nodes, u)
            end
        end
    end
    fill!(c.neigh_pos, 0)
    renumber_labels!(label, c.neigh_pos)
    label
end

"""Type to record neighbor labels and their counts."""
type NeighComm
  neigh_pos::Vector{Int}
  neigh_cnt::Vector{Float64}
  neigh_last::Int
end

"""Fast shuffle Array `a` in UnitRange `r` inplace."""
function range_shuffle!(r::UnitRange, a::AbstractVector)
    (r.start > 0 && r.stop <= length(a)) || error("out of bounds")
    @inbounds for i=length(r):-1:2
        j = rand(1:i)
        ii = i + r.start - 1
        jj = j + r.start - 1
        a[ii],a[jj] = a[jj],a[ii]
    end
end

"""Return the most frequency label."""
function vote!(g, weights, m, c, u, lambda)
    @inbounds for i=1:c.neigh_last-1
        c.neigh_cnt[c.neigh_pos[i]] = -1
    end
    c.neigh_last = 1
    c.neigh_pos[1] = m[u]
    c.neigh_cnt[c.neigh_pos[1]] = 0
    c.neigh_last = 2
    max_cnt = 0.
    for e in out_edges(u,g)
      w = weights[edge_index(e)]
      neigh = target(e)
      if u == neigh
        	w *= lambda/2
      end
        neigh_comm = m[neigh]
        if c.neigh_cnt[neigh_comm] < 0
            c.neigh_cnt[neigh_comm] = 0
            c.neigh_pos[c.neigh_last] = neigh_comm
            c.neigh_last += 1
        end
        c.neigh_cnt[neigh_comm] += w
        if c.neigh_cnt[neigh_comm] > max_cnt
          max_cnt = c.neigh_cnt[neigh_comm]
        end
    end
    # ties breaking randomly
    # range_shuffle!(1:c.neigh_last-1, c.neigh_pos)
    for lbl in c.neigh_pos
      if c.neigh_cnt[lbl] == max_cnt
        return lbl
      end
    end
end

function renumber_labels!(membership::Vector{Int}, label_counters::Vector{Int})
    N = length(membership)
    (maximum(membership) > N || minimum(membership) < 1) && error("Label must between 1 and |V|")
    j = 1
    @inbounds for i=1:length(membership)
        k = membership[i]
        if k >= 1
            if label_counters[k] == 0
                # We have seen this label for the first time
                label_counters[k] = j
                k = j
                j += 1
            else
                k = label_counters[k]
            end
        end
        membership[i] = k
    end
end

function collapse_graph{V}(g::AbstractGraph{V}, membership::Vector{Int},
    collapsed_edge_weights::Vector{Dict{Int,Float64}}, collapsed_weights::Vector{Float64})
    nb_comm = maximum(membership)

    resize!(collapsed_edge_weights, nb_comm)
    @inbounds for i=1:nb_comm
        collapsed_edge_weights[i] = Dict{Int,Float64}()
    end

    for e in edges(g)
        u = source(e,g)
        v = target(e,g)
        u_idx = vertex_index(u,g)
        v_idx = vertex_index(v,g)
        u_comm = membership[u_idx]
        v_comm = membership[v_idx]

        # for special case of undirected network
        if !is_directed(g)
            u_comm, v_comm = minmax(u_comm, v_comm)
        end

        if haskey(collapsed_edge_weights[u_comm], v_comm)
            collapsed_edge_weights[u_comm][v_comm] += 1
        else
            collapsed_edge_weights[u_comm][v_comm] = 1
        end
    end

    collapsed_graph = simple_graph(nb_comm, is_directed=false)
    empty!(collapsed_weights)

    @inbounds for u=1:nb_comm
        for (v,w) in collapsed_edge_weights[u]
            add_edge!(collapsed_graph, u, v)
            push!(collapsed_weights, w)
        end
    end

    collapsed_graph
end

function collapse_graph{V}(g::AbstractGraph{V}, membership::Vector{Int}, edge_weights::Vector{Float64},
    collapsed_edge_weights::Vector{Dict{Int,Float64}}, collapsed_weights::Vector{Float64})
    nb_comm = maximum(membership)

    resize!(collapsed_edge_weights, nb_comm)
    @inbounds for i=1:nb_comm
        collapsed_edge_weights[i] = Dict{Int,Float64}()
    end

    for e in edges(g)
    	w = edge_weights[edge_index(e,g)]
        u = source(e,g)
        v = target(e,g)
        u_idx = vertex_index(u,g)
        v_idx = vertex_index(v,g)
        u_comm = membership[u_idx]
        v_comm = membership[v_idx]

        # for special case of undirected network
        if !is_directed(g)
            u_comm, v_comm = minmax(u_comm, v_comm)
        end

        if haskey(collapsed_edge_weights[u_comm], v_comm)
            collapsed_edge_weights[u_comm][v_comm] += w
        else
            collapsed_edge_weights[u_comm][v_comm] = w
        end
    end

    collapsed_graph = simple_graph(nb_comm, is_directed=false)
    empty!(collapsed_weights)

    @inbounds for u=1:nb_comm
        for (v,w) in collapsed_edge_weights[u]
            add_edge!(collapsed_graph, u, v)
            push!(collapsed_weights, w)
        end
    end

    collapsed_graph
end

function from_coarser_partition!(membership::Vector{Int}, coarser_membership::Vector{Int})
    @inbounds for u=1:length(membership)
        # what is the community of the node
        u_comm_level1 = membership[u]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_membership[u_comm_level1]
        membership[u] = u_comm_level2
    end
end
