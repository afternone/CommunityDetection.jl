type DiFlowGroup{V}
    nodes::Set{V}
    inner_prob::Float64
    exit_prob::Float64
    iso_prob::Float64
end

type DiFlowPartition{V} <: AbstractPartition{V}
    flowgraph::DiFlowGraph{V}
    membership::Vector{Int}
    community::Dict{Int,DiFlowGroup{V}}
    total_exit_prob::Float64
end

# require interface
graph{V}(fp::DiFlowPartition{V}) = fp.flowgraph.graph
membership{V}(fp::DiFlowPartition{V}) = fp.membership
membership{V}(fp::DiFlowPartition{V}, u::V) = fp.membership[vertex_index(u, fp.flowgraph.graph)]
community{V}(fp::DiFlowPartition{V}) = keys(fp.community)
is_right_direction{V}(fp::DiFlowPartition{V}, quality, new_quality) = new_quality < quality
still_running{V}(fp::DiFlowPartition{V}, once_diff, ϵ) = once_diff < -ϵ

# construction
function diflow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(g)); τ=0.15)
    fg = diflow_graph(g, weights, τ=τ)
    diflow_partition(fg)
end

function diflow_partition{V}(fg::DiFlowGraph{V})
    n = num_vertices(fg.graph)
    fp = DiFlowPartition{V}(fg, collect(1:n), Dict{Int,DiFlowGroup{V}}(), 0.0)
    update_partition!(fp)
    fp
end

# mutation
function update_partition!{V}(partition::DiFlowPartition{V})
    maximum(partition.membership) ≤ num_vertices(partition.flowgraph.graph) || error("maximum(membership) must less than num_vertices(g)")
    minimum(partition.membership) > 0 || error("value of membership must be positive integer")

    fg = partition.flowgraph
    g = fg.graph

    n = num_vertices(g)
    empty!(partition.community)

    for u in vertices(g)
        u_idx = vertex_index(u, g)
        comm_idx = partition.membership[u_idx]
        if haskey(partition.community, comm_idx)
            ni = length(partition.community[comm_idx].nodes)
            partition.community[comm_idx].exit_prob += fg.tau*(n-ni-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[comm_idx].inner_prob/n -
                (1-fg.tau)*partition.community[comm_idx].iso_prob/n
            push!(partition.community[comm_idx].nodes, u)
            partition.community[comm_idx].inner_prob += fg.visit_prob[u_idx]
            if out_degree(u, g) > 0
                for e in out_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = target(e, g)
                    v_idx = vertex_index(v, g)
                    if !in(v, partition.community[comm_idx].nodes)
                        partition.community[comm_idx].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    end
                end
                for e in in_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = source(e, g)
                    v_idx = vertex_index(v, g)
                    if in(v, partition.community[comm_idx].nodes)
                        partition.community[comm_idx].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    end
                end
            else
                for e in in_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = source(e, g)
                    v_idx = vertex_index(v, g)
                    if in(v, partition.community[comm_idx].nodes)
                        partition.community[comm_idx].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    end
                end
                partition.community[comm_idx].iso_prob += fg.visit_prob[u_idx]
                partition.community[comm_idx].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-ni-1)/n
            end
        else
            if out_degree(u, fg.graph) > 0
                exit_prob = fg.tau*fg.visit_prob[u_idx]*(n-1)/n
                for e in out_edges(u, fg.graph)
                    e_idx = edge_index(e, fg.graph)
                    exit_prob += (1.0-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
                partition.community[comm_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob, 0.0)
            else
                partition.community[comm_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], fg.visit_prob[u_idx]*(n-1)/n, fg.visit_prob[u_idx])
            end
        end
    end

    partition.total_exit_prob = 0.0
    for group in values(partition.community)
        partition.total_exit_prob += group.exit_prob
    end
end

"Move a node to a new community and update the partition, this also removes any empty communities."
function move_node!{V}(partition::DiFlowPartition{V}, u::V, new_comm::Int)
    haskey(partition.community, new_comm) || error("partition has no community $new_comm")

    fg = partition.flowgraph
    g = fg.graph
    n = num_vertices(g)
    u_idx = vertex_index(u, g)
    old_comm = partition.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
        n_old = length(partition.community[old_comm].nodes)
        n_new = length(partition.community[new_comm].nodes)


        partition.community[old_comm].exit_prob += -fg.tau*(n-n_old+1)*fg.visit_prob[u_idx]/n +
                fg.tau*partition.community[old_comm].inner_prob/n +
                (1-fg.tau)*partition.community[old_comm].iso_prob/n
        partition.total_exit_prob += -fg.tau*(n-n_old+1)*fg.visit_prob[u_idx]/n +
                fg.tau*partition.community[old_comm].inner_prob/n +
                (1-fg.tau)*partition.community[old_comm].iso_prob/n
        partition.community[new_comm].exit_prob += fg.tau*(n-n_new-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[new_comm].inner_prob/n -
                (1-fg.tau)*partition.community[new_comm].iso_prob/n
        partition.total_exit_prob += fg.tau*(n-n_new-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[new_comm].inner_prob/n -
                (1-fg.tau)*partition.community[new_comm].iso_prob/n

        # remove node from old community
        delete!(partition.community[old_comm].nodes, u)
        # change of visit probability of the old community
        partition.community[old_comm].inner_prob -= partition.flowgraph.visit_prob[u_idx]

        # add node to new community
        push!(partition.community[new_comm].nodes, u)
        # change of inner probability of the new community
        partition.community[new_comm].inner_prob += partition.flowgraph.visit_prob[u_idx]

        if out_degree(u, g) > 0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                #if !in(v, partition.community[new_comm].nodes)
                if v_comm != new_comm
                    partition.community[new_comm].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
                #if !in(v, partition.community[old_comm].nodes)
                if v_comm != old_comm
                    partition.community[old_comm].exit_prob -= (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob -= (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
            end
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                #if in(v, partition.community[new_comm].nodes)
                if v_comm == new_comm
                    partition.community[new_comm].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
                #if in(v, partition.community[old_comm].nodes)
                if v_comm == old_comm
                    partition.community[old_comm].exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
            end
        else
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                #if in(v, partition.community[new_comm].nodes)
                if v_comm == new_comm
                    partition.community[new_comm].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
                #if in(v, partition.community[old_comm].nodes)
                if v_comm == old_comm
                    partition.community[old_comm].exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.taotal_exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
            end
            partition.community[new_comm].iso_prob += fg.visit_prob[u_idx]
            partition.community[old_comm].iso_prob -= fg.visit_prob[u_idx]
            partition.community[new_comm].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-n_new-1)/n
            partition.total_exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-n_new-1)/n
            partition.community[old_comm].exit_prob += -(1-fg.tau)*fg.visit_prob[u_idx]*(n-n_old+1)/n
            partition.total_exit_prob += -(1-fg.tau)*fg.visit_prob[u_idx]*(n-n_old+1)/n
        end

        # if the old community is empty after remove node u, we remove it
        if isempty(partition.community[old_comm].nodes)
            delete!(partition.community, old_comm)
        end

        # update the membership vector
        partition.membership[u_idx] = new_comm
    end
end

# retrun p*log(p)
plogp(x) = x > 0.0 ? x*log(x) : 0.0
plogp(xs::Vector{Float64}) = Float64[plogp(x) for x in xs]

"Returns the difference in average decribe length if we move a node to a new community"
function diff_move{V}(partition::DiFlowPartition{V}, u::V, new_comm::Int)
    u_idx = vertex_index(u, partition.flowgraph.graph)
    old_comm = partition.membership[u_idx]
    fg = partition.flowgraph
    g = fg.graph
    n = num_vertices(g)
    δL = 0.0

    #println(partition.community[old_comm].exit_prob)

    if new_comm != old_comm
        δexit_prob_old_comm = 0.0
        δexit_prob_new_comm = 0.0
        n_old = length(partition.community[old_comm].nodes)
        n_new = length(partition.community[new_comm].nodes)

        δexit_prob_old_comm += -fg.tau*(n-n_old+1)*fg.visit_prob[u_idx]/n +
                fg.tau*partition.community[old_comm].inner_prob/n +
                (1-fg.tau)*partition.community[old_comm].iso_prob/n
        δexit_prob_new_comm += fg.tau*(n-n_new-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[new_comm].inner_prob/n -
                (1-fg.tau)*partition.community[new_comm].iso_prob/n

        if out_degree(u, g) > 0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                #if !in(v, partition.community[new_comm].nodes)
                if v_comm != new_comm
                    δexit_prob_new_comm += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
                #if !in(v, partition.community[old_comm].nodes)
                if v_comm != old_comm
                    δexit_prob_old_comm -= (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
            end
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                #if in(v, partition.community[new_comm].nodes)
                if v_comm == new_comm
                    δexit_prob_new_comm -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
                #if in(v, partition.community[old_comm].nodes)
                if v_comm == old_comm
                    δexit_prob_old_comm += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
            end
        else
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                #if in(v, partition.community[new_comm].nodes)
                if v_comm == new_comm
                    δexit_prob_new_comm -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
                #if in(v, partition.community[old_comm].nodes)
                if v_comm == old_comm
                    δexit_prob_old_comm += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
            end
            δexit_prob_new_comm += (1-fg.tau)*fg.visit_prob[u_idx]*(n-n_new-1)/n
            δexit_prob_old_comm += -(1-fg.tau)*fg.visit_prob[u_idx]*(n-n_old+1)/n
        end

        δtotal_exit_prob = δexit_prob_old_comm + δexit_prob_new_comm

        δL1 = plogp(partition.total_exit_prob + δtotal_exit_prob) - plogp(partition.total_exit_prob)
        δL2 = -2(plogp(partition.community[old_comm].exit_prob + δexit_prob_old_comm) - plogp(partition.community[old_comm].exit_prob))
        δL3 = -2(plogp(partition.community[new_comm].exit_prob + δexit_prob_new_comm) - plogp(partition.community[new_comm].exit_prob))
        δL4 = plogp(partition.community[old_comm].exit_prob + δexit_prob_old_comm + partition.community[old_comm].inner_prob - fg.visit_prob[u_idx]) -
            plogp(partition.community[old_comm].exit_prob + partition.community[old_comm].inner_prob)
        δL5 = plogp(partition.community[new_comm].exit_prob + δexit_prob_new_comm + partition.community[new_comm].inner_prob + fg.visit_prob[u_idx]) -
            plogp(partition.community[new_comm].exit_prob + partition.community[new_comm].inner_prob)
        δL = δL1 + δL2 + δL3 + δL4 + δL5
    end

    δL
end

"Give the average decribe length of the partition."
function quality{V}(partition::DiFlowPartition{V})
    L1 = plogp(partition.total_exit_prob)
    L2 = -2sum(plogp([partition.community[i].exit_prob for i in keys(partition.community)]))
    L3 = -sum(plogp(partition.flowgraph.visit_prob))
    L4 = sum(plogp([partition.community[i].exit_prob+partition.community[i].inner_prob for i in keys(partition.community)]))

    L1 + L2 + L3 + L4
end

"""
Creates a graph with communities as node and links as weights between communities.

The weight of the edges in the new graph is simply the sum of the weight
of the edges between the communities. The size of a node in the new graph
is simply the size of the community in the old graph.
"""
function collapse_partition{V}(partition::DiFlowPartition{V})
    fg = partition.flowgraph
    g = fg.graph
    num_comm = length(partition.community)
    collapsed_trans_prob = Dict{Int,Float64}[]
    for i=1:num_comm
        push!(collapsed_trans_prob, Dict{Int,Float64}())
    end

    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_idx = vertex_index(u, g)
        v_idx = vertex_index(v, g)
        u_comm = partition.membership[u_idx]
        v_comm = partition.membership[v_idx]

        # we just skip self loop
        if u_comm != v_comm
            if haskey(collapsed_trans_prob[u_comm], v_comm)
                collapsed_trans_prob[u_comm][v_comm] += fg.trans_prob[e_idx]
            else
                collapsed_trans_prob[u_comm][v_comm] = fg.trans_prob[e_idx]
            end
        end
    end

    graph = simple_graph(num_comm)
    graph_trans_prob = Float64[]
    graph_visit_prob = Array(Float64, num_comm)

    for u_comm=1:num_comm
        graph_visit_prob[u_comm] = partition.community[u_comm].inner_prob
        for (v_comm, trans_prob) in collapsed_trans_prob[u_comm]
            add_edge!(graph, u_comm, v_comm)
            push!(graph_trans_prob, trans_prob)
        end
    end
    graph_trans_prob = trans_prob_directed(graph, graph_trans_prob)
    fg = DiFlowGraph(graph, fg.tau, graph_visit_prob, graph_trans_prob)
    diflow_partition(fg)
end
