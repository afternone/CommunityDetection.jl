using Graphs

function max_label{V}(g::AbstractGraph{V}, membership::Vector{Int}, 
    weights::Vector{Float64}, u::V, lammda::Float64, dominant_labels::Vector{Int}, label_cnt::Dict{Int,Float64})
    empty!(label_cnt)
    empty!(dominant_labels)
    max_cnt = 0.0;
    for e in out_edges(u,g)
        e_idx = edge_index(e,g)
        w = weights[e_idx]
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_lbl = membership[v_idx]
        if u == v
        	w *= lammda/2
        end
        label_cnt[v_lbl] = get(label_cnt, v_lbl, 0.0) + w
        if label_cnt[v_lbl] > max_cnt
        	max_cnt = label_cnt[v_lbl]
        end
    end
    for (lbl, wei) in label_cnt
    	if wei >= max_cnt
    		push!(dominant_labels, lbl)
    	end
    end
    dominant_labels[rand(1:length(dominant_labels))]
end

function move_node!{V}(g::AbstractGraph{V}, membership::Vector{Int}, 
    weights::Vector{Float64}, u::V, lammda::Float64, dominant_labels::Vector{Int}, label_cnt::Dict{Int,Float64})
    running = false
    if out_degree(u, g) > 0
        old_lbl = membership[vertex_index(u,g)]
        new_lbl = max_label(g, membership, weights, u, lammda, dominant_labels, label_cnt)
        if new_lbl != old_lbl
            membership[vertex_index(u,g)] = new_lbl
            running = true
        end
    end
    running
end

function hlpa{V}(g::AbstractGraph{V}, weights::Vector{Float64}, lammda::Float64, 
    dominant_labels::Vector{Int}, label_cnt::Dict{Int,Float64})
    label = collect(1:num_vertices(g))
    runing_nodes = Set(collect(vertices(g)))

    while !isempty(runing_nodes)
        order = shuffle(collect(runing_nodes))
        for u in order
            if !move_node!(g, label, weights, u, lammda, dominant_labels, label_cnt)
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

function hlpa{V}(g::AbstractGraph{V}; lammda::Float64=1.0)
    dominant_labels = Int[]
    sizehint!(dominant_labels, 100)
    label_cnt = Dict{Int, Float64}()
    sizehint!(label_cnt, 200)
    collapsed_edge_weights = Dict{Int,Float64}[]
    sizehint!(collapsed_edge_weights, num_vertices(g))
    collapsed_weights = Float64[]
    sizehint!(collapsed_weights, num_edges(g))
    s = similarity(g)
    m = hlpa(g, s, lammda, dominant_labels, label_cnt)
    nb_comm1 = maximum(m)
    h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
    m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m
end

function hlpa{V}(g::AbstractGraph{V}, s::Vector{Float64}; lammda::Float64=1.0)
    dominant_labels = Int[]
    sizehint!(dominant_labels, 100)
    label_cnt = Dict{Int, Float64}()
    sizehint!(label_cnt, 200)
    collapsed_edge_weights = Dict{Int,Float64}[]
    sizehint!(collapsed_edge_weights, num_vertices(g))
    collapsed_weights = Float64[]
    sizehint!(collapsed_weights, num_edges(g))
    m = hlpa(g, s, lammda, dominant_labels, label_cnt)
    nb_comm1 = maximum(m)
    h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
    m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m
end

function hlpa_record{V}(g::AbstractGraph{V}; lammda::Float64=1.0)
    label_rec = Vector{Int}[]
    graph_rec = typeof(g)[]
    ew_rec = Vector{Float64}[]

    dominant_labels = Int[]
    sizehint!(dominant_labels, 100)
    label_cnt = Dict{Int, Float64}()
    sizehint!(label_cnt, 200)
    collapsed_edge_weights = Dict{Int,Float64}[]
    sizehint!(collapsed_edge_weights, num_vertices(g))
    collapsed_weights = Float64[]
    sizehint!(collapsed_weights, num_edges(g))

    s = similarity(g)
    m = hlpa(g, s, lammda, dominant_labels, label_cnt)
    nb_comm1 = maximum(m)
    h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
    push!(label_rec, copy(m))
    push!(graph_rec, h)
    push!(ew_rec, copy(collapsed_weights))

    m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        push!(label_rec, m1)
        push!(graph_rec, h)
        push!(ew_rec, copy(collapsed_weights))

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m, label_rec, graph_rec, ew_rec
end

function hlpa_record{V}(g::AbstractGraph{V}, s::Vector{Float64}; lammda::Float64=1.0)
    label_rec = Vector{Int}[]
    graph_rec = typeof(g)[]
    ew_rec = Vector{Float64}[]

    dominant_labels = Int[]
    sizehint!(dominant_labels, 100)
    label_cnt = Dict{Int, Float64}()
    sizehint!(label_cnt, 200)
    collapsed_edge_weights = Dict{Int,Float64}[]
    sizehint!(collapsed_edge_weights, num_vertices(g))
    collapsed_weights = Float64[]
    sizehint!(collapsed_weights, num_edges(g))

    m = hlpa(g, s, lammda, dominant_labels, label_cnt)
    nb_comm1 = maximum(m)
    h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
    push!(label_rec, copy(m))
    push!(graph_rec, h)
    push!(ew_rec, copy(collapsed_weights))

    m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        push!(label_rec, m1)
        push!(graph_rec, h)
        push!(ew_rec, copy(collapsed_weights))

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h = collapse_graph(g, m, collapsed_edge_weights, collapsed_weights)
        

        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m, label_rec, graph_rec, ew_rec
end

function hlpa1{V}(g::AbstractGraph{V}; lammda::Float64=1.0)
    dominant_labels = Int[]
    sizehint!(dominant_labels, 100)
    label_cnt = Dict{Int, Float64}()
    sizehint!(label_cnt, 200)
    collapsed_edge_weights = Dict{Int,Float64}[]
    sizehint!(collapsed_edge_weights, num_vertices(g))
    collapsed_weights = Float64[]
    sizehint!(collapsed_weights, num_edges(g))
    s = similarity(g)
    m = hlpa(g, s, lammda, dominant_labels, label_cnt)
    nb_comm1 = maximum(m)
    h = collapse_graph(g, m, s, collapsed_edge_weights, collapsed_weights)
    m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
    from_coarser_partition!(m, m1)
    nb_comm2 = maximum(m)

    while nb_comm1 > nb_comm2
        nb_comm1 = nb_comm2
        h = collapse_graph(g, m, s, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb1 = maximum(m)
        h = collapse_graph(g, m, s, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb2 = maximum(m)
        h = collapse_graph(g, m, s, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb3 = maximum(m)
        h = collapse_graph(g, m, s, collapsed_edge_weights, collapsed_weights)
        m1 = hlpa(h, collapsed_weights, lammda, dominant_labels, label_cnt)
        from_coarser_partition!(m, m1)
        nb4 = maximum(m)
        nb_comm2 = (nb1+nb2+nb3+nb4)/4
    end
    m
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
