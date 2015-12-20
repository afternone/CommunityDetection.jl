function label_count{V}(g::AbstractGraph{V}, membership::Vector{Int}, cn::Vector{Int}, u::V)
    label_cnt = Dict{Int, Int}()
    for e in out_edges(u,g)
        e_idx = edge_index(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_lbl = membership[v_idx]
        if haskey(label_cnt, v_lbl)
            label_cnt[v_lbl] += 1 + cn[e_idx]
        else
            label_cnt[v_lbl] = 1 + cn[e_idx]
        end
    end
    label_cnt
end

function label_count{V}(g::AbstractGraph{V}, membership::Vector{Int}, u::V)
    label_cnt = Dict{Int, Int}()
    for v in out_neighbors(u,g)
        v_idx = vertex_index(v,g)
        v_lbl = membership[v_idx]
        if haskey(label_cnt, v_lbl)
            label_cnt[v_lbl] += 1
        else
            label_cnt[v_lbl] = 1
        end
    end
    label_cnt
end

function move_node!{V}(g::AbstractGraph{V}, membership::Vector{Int}, cn::Vector{Int}, u::V)
    old_lbl = membership[vertex_index(u,g)]
    label_cnt = label_count(g, membership, cn, u)
    neigh_lbls = collect(keys(label_cnt))
    neigh_lbls_cnt = collect(values(label_cnt))
    order = collect(1:length(label_cnt))
    shuffle!(order)
    max_cnt = 0
    max_lbl = old_lbl
    running = false
    for i in order
        if neigh_lbls_cnt[i] > max_cnt
            max_cnt = neigh_lbls_cnt[i]
            max_lbl = neigh_lbls[i]
        end
    end
    if max_lbl != old_lbl
        membership[vertex_index(u,g)] = max_lbl
        running = true
    end
    running
end

function move_node!{V}(g::AbstractGraph{V}, membership::Vector{Int}, u::V)
    old_lbl = membership[vertex_index(u,g)]
    label_cnt = label_count(g, membership, u)
    neigh_lbls = collect(keys(label_cnt))
    neigh_lbls_cnt = collect(values(label_cnt))
    order = collect(1:length(label_cnt))
    shuffle!(order)
    max_cnt = 0
    max_lbl = old_lbl
    running = false
    for i in order
        if neigh_lbls_cnt[i] > max_cnt
            max_cnt = neigh_lbls_cnt[i]
            max_lbl = neigh_lbls[i]
        end
    end
    if max_lbl != old_lbl
        membership[vertex_index(u,g)] = max_lbl
        running = true
    end
    running
end

function nsdlpa{V}(g::AbstractGraph{V})
    label = collect(1:num_vertices(g))
    runing_nodes = Set(collect(vertices(g)))
    cn = num_common_neighbors(g)

    while !isempty(runing_nodes)
        order = shuffle(collect(runing_nodes))
        for u in order
            if !move_node!(g, label, cn, u)
                delete!(runing_nodes, u)
            else
                for v in out_neighbors(u,g)
                    push!(runing_nodes, v)
                end
            end
        end
    end
    label
end

function lpa{V}(g::AbstractGraph{V})
    label = collect(1:num_vertices(g))
    runing_nodes = Set(collect(vertices(g)))

    while !isempty(runing_nodes)
        order = shuffle(collect(runing_nodes))
        for u in order
            if !move_node!(g, label, u)
                delete!(runing_nodes, u)
            else
                for v in out_neighbors(u,g)
                    push!(runing_nodes, v)
                end
            end
        end
    end
    label
end
