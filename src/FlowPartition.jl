# type of group
type FlowGroup{V}
    nodes::Set{V}
    inner_prob::Float64
    exit_prob::Float64
end

# type of flow partition
type FlowPartition{V} <: AbstractPartition{V}
    flowgraph::FlowGraph{V}
    membership::Vector{Int}
    community::Dict{Int,FlowGroup{V}}
    total_exit_prob::Float64
end

# require interface
graph{V}(fp::FlowPartition{V}) = fp.flowgraph.graph
membership{V}(fp::FlowPartition{V}) = fp.membership
membership{V}(fp::FlowPartition{V}, u::V) = fp.membership[vertex_index(u, fp.flowgraph.graph)]
community{V}(fp::FlowPartition{V}) = keys(fp.community)
is_right_direction{V}(fp::FlowPartition{V}, quality, new_quality) = new_quality < quality
still_running{V}(fp::FlowPartition{V}, once_diff, ϵ) = once_diff < -ϵ

# construction
function flow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(g)))
    fg = flow_graph(g, weights)
    flow_partition(fg)
end

function flow_partition{V}(fg::FlowGraph{V})
    n = num_vertices(fg.graph)
    fp = FlowPartition{V}(fg, collect(1:n), Dict{Int,FlowGroup{V}}(), 0.0)
    update_partition!(fp)
    fp
end

# mutation
"update partition when membership vector changed"
function update_partition!{V}(fp::FlowPartition{V})
    membership = copy(fp.membership)
    update_partition!(fp, membership)
end

"update partition based on membership vector"
function update_partition!{V}(fp::FlowPartition{V})
    maximum(fp.membership) ≤ num_vertices(fp.flowgraph.graph) || error("maximum(membership) must less than num_vertices(graph)")
    minimum(fp.membership) > 0 || error("value of membership must be positive integer")

    fg = fp.flowgraph
    g = fg.graph
    # update membership vector and communities
    empty!(fp.community)
    #fp.community = Dict{Int,FlowGroup{V}}()
    for u in vertices(g)
        u_idx = vertex_index(u, g)
        comm_idx = fp.membership[u_idx]
        if haskey(fp.community, comm_idx)
            push!(fp.community[comm_idx].nodes, u)
            fp.community[comm_idx].inner_prob += fg.visit_prob[u_idx]
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                if !in(v, fp.community[comm_idx].nodes)
                    fp.community[comm_idx].exit_prob += fg.trans_prob[e_idx]
                else
                    fp.community[comm_idx].exit_prob -= fg.trans_prob[e_idx]
                end
            end
        else
            exit_prob = 0.0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                exit_prob += fg.trans_prob[e_idx]
            end
            fp.community[comm_idx] = FlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob)
        end
    end

    fp.total_exit_prob = 0.0
    for group in values(fp.community)
        fp.total_exit_prob += group.exit_prob
    end
end

"Move a node to a new community and update the partition, this also removes any empty communities."
function move_node!{V}(fp::FlowPartition{V}, u::V, new_comm::Int)
    haskey(fp.community, new_comm) || error("partition has no community $new_comm")
    fg = fp.flowgraph
    g = fg.graph

    u_idx = vertex_index(u, g)
    old_comm = fp.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
        # remove node from old community
        delete!(fp.community[old_comm].nodes, u)
        # change of visit probability of the old community
        fp.community[old_comm].inner_prob -= fg.visit_prob[u_idx]

        # add node to new community
        push!(fp.community[new_comm].nodes, u)
        # change of inner probability of the new community
        fp.community[new_comm].inner_prob += fg.visit_prob[u_idx]

        for e in out_edges(u, g)
            e_idx = edge_index(e, g)
            v = target(e, g)
            v_idx = vertex_index(v,g)
            v_comm = fp.membership[v_idx]
            # change of exit probability of the old community
            #if !in(v, partition.community[old_comm].nodes)
            if v_comm != old_comm
                fp.community[old_comm].exit_prob -= fg.trans_prob[e_idx]
                fp.total_exit_prob -= fg.trans_prob[e_idx]
            else
                fp.community[old_comm].exit_prob += fg.trans_prob[e_idx]
                fp.total_exit_prob += fg.trans_prob[e_idx]
            end
            # change of exit probability of the new community
            #if !in(v, partition.community[new_comm].nodes)
            if v_comm != new_comm
                fp.community[new_comm].exit_prob += fg.trans_prob[e_idx]
                fp.total_exit_prob += fg.trans_prob[e_idx]
            else
                fp.community[new_comm].exit_prob -= fg.trans_prob[e_idx]
                fp.total_exit_prob -= fg.trans_prob[e_idx]
            end
        end

        # if the old community is empty after remove node u, we remove it
        if isempty(fp.community[old_comm].nodes)
            delete!(fp.community, old_comm)
        end

        # update the membership vector
        fp.membership[u_idx] = new_comm
    end
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
            e_idx = edge_index(e, g)
            v = target(e, g)
            v_idx = vertex_index(v,g)
            v_comm = fp.membership[v_idx]
            #if in(v, fp.community[old_comm].nodes)
            if v_comm == old_comm
                δexit_prob_old_comm += fg.trans_prob[e_idx]
            else
                δexit_prob_old_comm -= fg.trans_prob[e_idx]
            end
            #if in(v, fp.community[new_comm].nodes)
            if v_comm == new_comm
                δexit_prob_new_comm -= fg.trans_prob[e_idx]
            else
                δexit_prob_new_comm += fg.trans_prob[e_idx]
            end
        end
        δtotal_exit_prob = δexit_prob_old_comm + δexit_prob_new_comm

        δL1 = plogp(fp.total_exit_prob + δtotal_exit_prob) - plogp(fp.total_exit_prob)
        δL2 = -2(plogp(fp.community[old_comm].exit_prob + δexit_prob_old_comm) - plogp(fp.community[old_comm].exit_prob))
        δL3 = -2(plogp(fp.community[new_comm].exit_prob + δexit_prob_new_comm) - plogp(fp.community[new_comm].exit_prob))
        δL4 = plogp(fp.community[old_comm].exit_prob + δexit_prob_old_comm + fp.community[old_comm].inner_prob - fg.visit_prob[u_idx]) -
            plogp(fp.community[old_comm].exit_prob + fp.community[old_comm].inner_prob)
        δL5 = plogp(fp.community[new_comm].exit_prob + δexit_prob_new_comm + fp.community[new_comm].inner_prob + fg.visit_prob[u_idx]) -
            plogp(fp.community[new_comm].exit_prob + fp.community[new_comm].inner_prob)
        δL = δL1 + δL2 + δL3 + δL4 + δL5
    end

    δL
end

"Give the average decribe length of the partition."
function quality{V}(fp::FlowPartition{V})
    L1 = plogp(fp.total_exit_prob)
    L2 = -2sum(plogp([fp.community[i].exit_prob for i in keys(fp.community)]))
    L3 = -sum(plogp(fp.flowgraph.visit_prob))
    L4 = sum(plogp([fp.community[i].exit_prob+fp.community[i].inner_prob for i in keys(fp.community)]))

    L1 + L2 + L3 + L4
end

"""
Creates a graph with communities as node and links as weights between communities.

The weight of the edges in the new graph is simply the sum of the weight
of the edges between the communities. The size of a node in the new graph
is simply the size of the community in the old graph.
"""
function collapse_partition{V}(fp::FlowPartition{V})
    fg = fp.flowgraph
    g = fg.graph
    num_comm = length(fp.community)

    collapsed_trans_prob = Array(Dict{Int,Float64}, num_comm)
    for i=1:num_comm
        collapsed_trans_prob[i] = Dict{Int,Float64}()
    end

    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_idx = vertex_index(u, g)
        v_idx = vertex_index(v, g)
        u_comm = fp.membership[u_idx]
        v_comm = fp.membership[v_idx]

        # because network is undirected, we always set u_comm < v_comm
        u_comm, v_comm = minmax(u_comm, v_comm)
        # we just skip self loop
        if u_comm != v_comm
            if haskey(collapsed_trans_prob[u_comm], v_comm)
                collapsed_trans_prob[u_comm][v_comm] += fg.trans_prob[e_idx]
            else
                collapsed_trans_prob[u_comm][v_comm] = fg.trans_prob[e_idx]
            end
        end
    end

    graph = simple_graph(num_comm, is_directed=false)
    graph_trans_prob = Float64[]
    graph_visit_prob = Array(Float64, num_comm)

    for u_comm=1:num_comm
        graph_visit_prob[u_comm] = fp.community[u_comm].inner_prob
        for (v_comm, trans_prob) in collapsed_trans_prob[u_comm]
            add_edge!(graph, u_comm, v_comm)
            push!(graph_trans_prob, trans_prob)
        end
    end

    fg = FlowGraph(graph, graph_visit_prob, graph_trans_prob)
    flow_partition(fg)
end
