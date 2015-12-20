using Graphs

# type of undirected flow graph
type FlowGraph{V}
    graph::AbstractGraph{V}
    visit_prob::Vector{Float64} # nodes visit probability
    edge_weights::AbstractEdgePropertyInspector{Float64} # relative edge weights
end

function flow_graph{V}(g::AbstractGraph{V})
    is_directed(g) && error("graph must undirected")
    m = num_edges(g)
    visit_prob = Float64[out_degree(u,g)/2m for u in vertices(g)]
    edge_weights = ConstantEdgePropertyInspector(1/2m)
    FlowGraph(g, visit_prob, edge_weights)
end

function flow_graph{V,T}(g::AbstractGraph{V}, weights::Vector{T})
    !is_directed(g) || error("graph must undirected")
    length(weights) == num_edges(g) || error("wrong edge weights length")

    total_weight = sum(weights)
    edge_weights = VectorEdgePropertyInspector(weights/2total_weight)
    visit_prob = zeros(num_vertices(g))
    for e in edges(g)
        u = source(e,g)
        v = target(e,g)
        w = edge_property(edge_weights, e, g)
        visit_prob[vertex_index(u,g)] += w
        visit_prob[vertex_index(v,g)] += w
    end
    FlowGraph(g, visit_prob, edge_weights)
end

function collapse_graph{V}(fg::FlowGraph{V}, membership::Vector{Int})
    length(membership) == num_vertices(fg.graph) || error("wrong membership length")
    g = fg.graph
    nb_comm = maximum(membership)
    collapsed_visit_prob = zeros(nb_comm)
    for u in vertices(g)
        u_idx = vertex_index(u,g)
        collapsed_visit_prob[membership[u_idx]] += fg.visit_prob[u_idx]
    end

    collapsed_edge_weights = [Dict{Int,Float64}() for i=1:nb_comm]
    for e in edges(g)
        u = source(e,g)
        v = target(e,g)
        u_comm = membership[vertex_index(u,g)]
        v_comm = membership[vertex_index(v,g)]
        # for unidrected
        u_comm, v_comm = minmax(u_comm, v_comm)
        # excluding self-loop
        if u_comm != v_comm
            collapsed_edge_weights[u_comm][v_comm] = get(collapsed_edge_weights[u_comm], v_comm, 0.0) + edge_property(fg.edge_weights, e, g)
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

    FlowGraph(collapsed_graph, collapsed_visit_prob, VectorEdgePropertyInspector(collapsed_weights))
end


# type of flow partition
type FlowPartition{V}
    flowgraph::FlowGraph{V}
    membership::Vector{Int}
    inner_prob::Vector{Float64}
    exit_prob::Vector{Float64}
    total_exit_prob::Float64
end

# construction
function flow_partition{V}(g::AbstractGraph{V})
    n = num_vertices(g)
    fp = FlowPartition{V}(flow_graph(g), collect(1:n), Float64[], Float64[], 0.0)
    update_partition!(fp)
    fp
end

function flow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T})
    n = num_vertices(g)
    fp = FlowPartition{V}(flow_graph(g, weights), collect(1:n), Float64[], Float64[], 0.0)
    update_partition!(fp)
    fp
end

function collapse_partition{V}(fp::FlowPartition{V})
    collapsed_fg = collapse_graph(fp.flowgraph, fp.membership)
    n = num_vertices(collapsed_fg.graph)
    collapsed_fp = FlowPartition(collapsed_fg, collect(1:n), Float64[], Float64[], 0.0)
    update_partition!(collapsed_fp)
    collapsed_fp
end

"update partition based on membership vector"
function update_partition!{V}(fp::FlowPartition{V})
    nb_comm = maximum(fp.membership)
    fg = fp.flowgraph
    g = fg.graph
    # reset partition
    fp.total_exit_prob = 0.0
    resize!(fp.inner_prob, nb_comm)
    resize!(fp.exit_prob, nb_comm)
    @inbounds for i=1:nb_comm
        fp.inner_prob[i] = 0.0
        fp.exit_prob[i] = 0.0
    end

    for u in vertices(g)
        u_idx = vertex_index(u,g)
        u_comm = fp.membership[u_idx]
        fp.inner_prob[u_comm] += fg.visit_prob[u_idx]
    end

    for e in edges(g)
        u = source(e,g)
        v = target(e,g)
        u_comm = fp.membership[vertex_index(u,g)]
        v_comm = fp.membership[vertex_index(v,g)]
        w = edge_property(fg.edge_weights, e, g)
        if v_comm != u_comm
            fp.exit_prob[u_comm] += w
            fp.exit_prob[v_comm] += w
            fp.total_exit_prob += 2w
        end
    end
end

"Move a node to a new community and update the partition, this also removes any empty communities."
function move_node!{V}(fp::FlowPartition{V}, u::V, new_comm::Int)
    fg = fp.flowgraph
    g = fg.graph

    u_idx = vertex_index(u, g)
    old_comm = fp.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
        # change of visit probability of the old community
        fp.inner_prob[old_comm] -= fg.visit_prob[u_idx]

        # change of inner probability of the new community
        fp.inner_prob[new_comm] += fg.visit_prob[u_idx]

        for e in out_edges(u, g)
            v = target(e, g)
            v_idx = vertex_index(v,g)
            v_comm = fp.membership[v_idx]
            w = edge_property(fg.edge_weights, e, g)
            # change of exit probability of the old community
            if v_comm != old_comm
                fp.exit_prob[old_comm] -= w
                fp.total_exit_prob -= w
            else
                fp.exit_prob[old_comm] += w
                fp.total_exit_prob += w
            end
            # change of exit probability of the new community
            #if !in(v, partition.community[new_comm].nodes)
            if v_comm != new_comm
                fp.exit_prob[new_comm] += w
                fp.total_exit_prob += w
            else
                fp.exit_prob[new_comm] -= w
                fp.total_exit_prob -= w
            end
        end

        # update the membership vector
        fp.membership[u_idx] = new_comm
    end
end

function renumber_communities!{V}(fp::FlowPartition{V})
    idx_map = Dict{Int,Int}()
    j = 1
    for i in fp.membership
        if !haskey(idx_map, i)
            idx_map[i] = j
            j += 1
        end
    end

    nb_comm = length(idx_map)
    new_inner_prob = Array(Float64, nb_comm)
    new_exit_prob = Array(Float64, nb_comm)
    for (old_comm_idx, new_comm_idx) in idx_map
        new_inner_prob[new_comm_idx] = fp.inner_prob[old_comm_idx]
        new_exit_prob[new_comm_idx] = fp.exit_prob[old_comm_idx]
    end

    @inbounds for i=1:num_vertices(fp.flowgraph.graph)
        fp.membership[i] = idx_map[fp.membership[i]]
    end
    fp.inner_prob = new_inner_prob
    fp.exit_prob = new_exit_prob
end

function from_coarser_partition!{V}(fp::FlowPartition{V}, coarser_fp::FlowPartition{V})
    for u=1:length(fp.membership)
        # what is the community of the node
        u_comm_level1 = fp.membership[u]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_fp.membership[u_comm_level1]
        fp.membership[u] = u_comm_level2
    end
    update_partition!(fp)
end

# retrun p*log(p)
plogp(x) = x > 0.0 ? x*log(x) : 0.0
plogp(xs::Vector{Float64}) = Float64[plogp(x) for x in xs]

"Returns the difference in average decribe length if we move a node to a new community"
function diff_move{V}(fp::FlowPartition{V}, u::V, new_comm::Int)
    fg = fp.flowgraph
    g = fg.graph
    u_idx = vertex_index(u, g)
    old_comm = fp.membership[u_idx]

    δL = 0.0
    if new_comm != old_comm
        δexit_prob_old_comm = 0.0
        δexit_prob_new_comm = 0.0
        for e in out_edges(u, g)
            w = edge_property(fg.edge_weights, e, g)
            v = target(e, g)
            v_idx = vertex_index(v,g)
            v_comm = fp.membership[v_idx]
            if v_comm == old_comm
                δexit_prob_old_comm += w
            else
                δexit_prob_old_comm -= w
            end
            if v_comm == new_comm
                δexit_prob_new_comm -= w
            else
                δexit_prob_new_comm += w
            end
        end
        δtotal_exit_prob = δexit_prob_old_comm + δexit_prob_new_comm

        δL1 = plogp(fp.total_exit_prob + δtotal_exit_prob) - plogp(fp.total_exit_prob)
        δL2 = -2(plogp(fp.exit_prob[old_comm] + δexit_prob_old_comm) - plogp(fp.exit_prob[old_comm]))
        δL3 = -2(plogp(fp.exit_prob[new_comm] + δexit_prob_new_comm) - plogp(fp.exit_prob[new_comm]))
        δL4 = plogp(fp.exit_prob[old_comm] + δexit_prob_old_comm + fp.inner_prob[old_comm] - fg.visit_prob[u_idx]) -
            plogp(fp.exit_prob[old_comm] + fp.inner_prob[old_comm])
        δL5 = plogp(fp.exit_prob[new_comm] + δexit_prob_new_comm + fp.inner_prob[new_comm] + fg.visit_prob[u_idx]) -
            plogp(fp.exit_prob[new_comm] + fp.inner_prob[new_comm])
        δL = δL1 + δL2 + δL3 + δL4 + δL5
    end

    δL
end

"Give the average decribe length of the partition."
function quality{V}(fp::FlowPartition{V})
    L1 = plogp(fp.total_exit_prob)
    L2 = -2sum(plogp(fp.exit_prob))
    L3 = -sum(plogp(fp.flowgraph.visit_prob))
    L4 = sum(plogp(fp.exit_prob+fp.inner_prob))

    L1 + L2 + L3 + L4
end

"get neighbor communities of node u"
function get_neigh_comms{V}(fp::FlowPartition{V}, u::V)
    g = fp.flowgraph.graph
    neigh_comms = Set{Int}()
    for v in out_neighbors(u, g)
        v_idx = vertex_index(v, g)
        push!(neigh_comms, fp.membership[v_idx])
    end
    neigh_comms
end

"Move nodes to other communities depending on how other communities are
considered."
function move_nodes!{V}(fp::FlowPartition{V}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    g = fp.flowgraph.graph
    memb = fp.membership
    itr = 0 # number of iterations
    total_diff = 0.0 # total difference while moving nodes
    once_diff = 2ϵ # difference for one loop
    n = num_vertices(g)
    num_moves = 2n # number of moved nodes during one loop
    vertex_order = collect(vertices(g))
    # As long as we keep on improving and we don't exceed the
    # maximum number of iterations and number of moves.
    while once_diff > ϵ && num_moves > n*δ && itr < max_itr
        # increase number of iterations
        itr += 1

        # initialize number of moves and imporvment
        num_moves = 0
        once_diff = 0.0
        shuffle!(vertex_order)

        # for each node
        for u in vertex_order
            # only take into account nodes of degree higher than zero
            if out_degree(u, g) > 0
                u_idx = vertex_index(u, g)
                # What is the current community of the node
                u_comm = memb[u_idx]
                # What is the difference per community if we move the node to one of
                # the other communities, and what is the extreme difference?
                min_diff = 0.0
                new_comm = u_comm


                # In which communities are its neighbours
                #neigh_comms = collect(get_neigh_comms(fp, u))
                neigh_comms = get_neigh_comms(fp, u)
                # tie break randomly, it has been found that this strategy can improv result
                #shuffle!(neigh_comms)
                # Loop through the communities of the neighbours
                for comm in neigh_comms
                    # Calculate the possible difference of the moving the node to that community
                    possible_diff = diff_move(fp, u, comm)
                    # We're only inserested in the minimum
                    if possible_diff < min_diff
                        min_diff = possible_diff
                        new_comm = comm
                    end
                end

                if new_comm != u_comm
                    # keep track of quality change
                    once_diff += min_diff
                    # actually move the node
                    move_node!(fp, u, new_comm)
                    num_moves += 1
                end
            end # if out_degree
        end # for

        # keep track of total difference over multiple loops
        total_diff += once_diff
    end # while
    renumber_communities!(fp)
    total_diff
end

"Optimize the provided partition."
function optimize_partition!{V}(fp::FlowPartition{V}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    # do one iteration of optimisation
    once_diff = move_nodes!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    # as long as there remains imprvoment iterate
    while once_diff < -ϵ
        # first create collapsed partition
        collapsed_fp = collapse_partition(fp)
        # Optimise partition for collapsed graph
        once_diff = move_nodes!(collapsed_fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
        # Make sure improvement on coarser scale is reflected on the
        # scale of the graph as a whole.
        from_coarser_partition!(fp, collapsed_fp)
    end

    # We renumber the communities to make sure we stick in the range
    # 1,...,r for r communities.
    # By default, we number the communities in decreasing order of size,
    # so that 1 is the largest community, 2 the second largest, etc...
    renumber_communities!(fp)
    # Return the quality of the current partition.
    quality(fp)
end

function find_partition!{V}(fp::FlowPartition{V}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    quality1 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    quality2 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)

    while quality2-quality1 < -ϵ
        quality1 = quality2
        quality2 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    end
    quality2
end

function find_partition!{V}(g::AbstractGraph{V}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    fp = flow_partition(g)
    quality1 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    quality2 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)

    while quality2-quality1 < -ϵ
        quality1 = quality2
        quality2 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    end
    fp
end

function find_partition!{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    fp = flow_partition(g, weights)
    quality1 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    quality2 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)

    while quality2-quality1 < -ϵ
        quality1 = quality2
        quality2 = optimize_partition!(fp, ϵ=ϵ, δ=δ, max_itr=max_itr)
    end
    fp
end
