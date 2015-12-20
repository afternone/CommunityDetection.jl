using Graphs

function label_count{V}(g::AbstractGraph{V}, membership::Vector{Int}, weights::Vector{Float64}, u::V)
    label_cnt = Dict{Int, Float64}()
    for e in unique(out_edges(u,g))
        e_idx = edge_index(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_lbl = membership[v_idx]
        label_cnt[v_lbl] = get(label_cnt, v_lbl, 0.0) + weights[e_idx]
    end
    label_cnt
end

function move_node!{V}(g::AbstractGraph{V}, membership::Vector{Int}, weights::Vector{Float64}, u::V)
    running = false
    if out_degree(u, g) > 0
        old_lbl = membership[vertex_index(u,g)]
        label_cnt = label_count(g, membership, weights, u)
        neigh_lbls = collect(keys(label_cnt))
        neigh_lbls_cnt = collect(values(label_cnt))
        order = collect(1:length(label_cnt))
        shuffle!(order)
        max_cnt = neigh_lbls_cnt[order[1]]
        max_lbl = neigh_lbls[order[1]]
        for i=2:length(order)
            if neigh_lbls_cnt[order[i]] > max_cnt
                max_cnt = neigh_lbls_cnt[order[i]]
                max_lbl = neigh_lbls[order[i]]
            end
        end
        if max_lbl != old_lbl
            membership[vertex_index(u,g)] = max_lbl
            running = true
        end
    end
    running
end

function permute_labels!(membership::Vector{Int})
    N = length(membership)
    if maximum(membership) > N || minimum(membership) < 1
        error("Label must between 1 and |V|")
    end
    label_counters = zeros(Int, N)
    j = 1
    for i=1:length(membership)
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

function lpa{V}(g::AbstractGraph{V}, weights::Vector{Float64})
    label = collect(1:num_vertices(g))
    runing_nodes = Set(collect(vertices(g)))

    while !isempty(runing_nodes)
        order = shuffle(collect(runing_nodes))
        for u in order
            if !move_node!(g, label, weights, u)
                delete!(runing_nodes, u)
            else
                for v in out_neighbors(u,g)
                    push!(runing_nodes, v)
                end
            end
        end
    end
    permute_labels!(label)
    label
end

function similarity{V,E}(g::AbstractGraph{V,E}, e::E)
    u = source(e, g)
    v = target(e, g)
    u_nei = out_neighbors(u, g)
    v_nei = out_neighbors(v, g)
    num_common_nei = length(intersect(u_nei, v_nei))
    # exclude u and v, so we minus 2
    num_total_nei = length(union(u_nei, v_nei)) - 2
    (num_common_nei + 1) / (num_total_nei + 1)
end

similarity{V,E}(g::AbstractGraph{V,E}) = Float64[similarity(g,e) for e in edges(g)]

function lpa{V}(g::AbstractGraph{V})
    s = similarity(g)
    m = lpa(g, s)
    nb_comm1 = maximum(m)
    h, w = collapse_graph(g, m)
    m1 = lpa(h, w)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h, w = collapse_graph(g, m)
        m1 = lpa(h, w)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h, w = collapse_graph(g, m)
        m1 = lpa(h, w)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h, w = collapse_graph(g, m)
        m1 = lpa(h, w)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h, w = collapse_graph(g, m)
        m1 = lpa(h, w)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m
end

function collapse_graph{V}(g::AbstractGraph{V}, membership::Vector{Int})
    nb_comm = maximum(membership)

    collapsed_edge_weights = Array(Dict{Int,Float64}, nb_comm)
    for i=1:nb_comm
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
    collapsed_weights = Float64[]

    for u=1:nb_comm
        for (v,w) in collapsed_edge_weights[u]
            add_edge!(collapsed_graph, u, v)
            push!(collapsed_weights, w)
        end
    end

    collapsed_graph, collapsed_weights
end

function from_coarser_partition!(membership::Vector{Int}, coarser_membership::Vector{Int})
    for u=1:length(membership)
        # what is the community of the node
        u_comm_level1 = membership[u]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_membership[u_comm_level1]
        membership[u] = u_comm_level2
    end
end
