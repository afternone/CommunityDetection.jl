function multi_greedy_merge!{V}(fp::FlowPartition{V}, ϵ::Float64=1e-5)
    fg = fp.flowgraph
    g = fg.graph
    edges_to_merge = edge_type(g)[]
    edges_diff = Float64[]

    for e in edges(g)
        u = source(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_comm = fp.membership[v_idx]
        _diff = diff_move(fp, u, v_comm)
        if _diff < -ϵ
            push!(edges_to_merge, e)
            push!(edges_diff, _diff)
        end
    end

    merged_nodes = Set{V}()
    perm_idx = sortperm(edges_diff)
    for i in perm_idx
        e = edges_to_merge[i]
        u = source(e,g)
        v = target(e,g)
        if !in(u, merged_nodes) && !in(v, merged_nodes)
            v_idx = vertex_index(v,g)
            v_comm = fp.membership[v_idx]
            move_node!(fp, u, v_comm)
            push!(merged_nodes, u)
            push!(merged_nodes, v)
        end
    end
    !isempty(edges_diff)
end

function multi_greedy!{V}(partition::AbstractPartition{V}, ϵ=1e-5; kargs...)
    # do one iteration of optimisation
    once_diff = move_nodes!(partition, ϵ; kargs...)
    collapsed_partition = collapse_partition(partition)
    running = multi_greedy_merge!(collapsed_partition, ϵ)

    while running || once_diff < -ϵ
        collapsed_partition = collapse_partition(partition)
        running = multi_greedy_merge!(collapsed_partition, ϵ)
        from_coarser_partition!(partition, collapsed_partition)
        once_diff = move_nodes!(partition, ϵ; kargs...)
    end

    renumber_communities!(partition)

    quality(partition)
end
