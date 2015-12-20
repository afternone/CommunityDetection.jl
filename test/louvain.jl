<<<<<<< HEAD
using Graphs

type MGraph{V}
    graph::AbstractGraph{V}
    strengths::Vector{Float64}
    self_weights::Vector{Float64}
=======
using Graphs,GraphPlot

# type of undirected modularity graph
type MGraph{V}
    graph::AbstractGraph{V}
    node_strengths_in::Vector{Float64}
    node_strengths_out::Vector{Float64}
    node_self_weights::Vector{Float64}
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
    edge_weights::AbstractEdgePropertyInspector{Float64}
    total_weight::Float64
end

function mgraph{V}(g::AbstractGraph{V})
<<<<<<< HEAD
    is_directed(g) && error("graph must undirected")
    m = num_edges(g)
    strengths = Float64[out_degree(u,g) for u in vertices(g)]
    edge_weights = ConstantEdgePropertyInspector(1.0)
    MGraph(g, strengths, zeros(num_vertices(g)), edge_weights, float(m))
end

function mgraph{V,T}(g::AbstractGraph{V}, weights::Vector{T})
    !is_directed(g) || error("graph must undirected")
    length(weights) == num_edges(g) || error("wrong edge weights length")

    total_weight = sum(weights)
    edge_weights = VectorEdgePropertyInspector(weights/2total_weight)
    strengths = zeros(num_vertices(g))
=======
    node_strengths_in = is_directed(g) ? Float64[in_degree(u,g) for u in vertices(g)] : Float64[]
    node_strengths_out = Float64[out_degree(u,g) for u in vertices(g)]
    edge_weights = ConstantEdgePropertyInspector(1.0)
    MGraph(g, node_strengths_in, node_strengths_out, zeros(num_vertices(g)), edge_weights, float(num_edges(g)))
end

function mgraph{V,T}(g::AbstractGraph{V}, weights::Vector{T})
    length(weights) == num_edges(g) || error("wrong edge weights length")

    total_weight = sum(weights)
    edge_weights = VectorEdgePropertyInspector(weights)
    node_strengths_in = is_directed(g) ? zeros(num_vertices(g)) : Float64[]
    node_strengths_out = zeros(num_vertices(g))

>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
    for e in edges(g)
        u = source(e,g)
        v = target(e,g)
        w = edge_property(edge_weights, e, g)
<<<<<<< HEAD
        strengths[vertex_index(u,g)] += w
        strengths[vertex_index(v,g)] += w
    end
    MGraph(g, strengths, zeros(num_vertices(g)), edge_weights, total_weight)
end

function collapse_graph{V}(mg::MGraph{V}, membership::Vector{Int})
    g = mg.graph
    length(membership) == num_vertices(g) || error("wrong membership length")

    nb_comm = maximum(membership)
    collapsed_strengths = zeros(nb_comm)
    collapsed_self_weights = zeros(nb_comm)
    for u in vertices(g)
        u_idx = vertex_index(u,g)
        collapsed_self_weights[membership[u_idx]] += mg.self_weights[u_idx]
    end

    collapsed_edge_weights = [Dict{Int,Float64}() for i=1:nb_comm]
    for e in edges(g)
        w = edge_property(mg.edge_weights, e, g)
=======
        node_strengths_out[vertex_index(u,g)] += w
        if is_directed(g)
            node_strengths_in[vertex_index(v,g)] += w
        else
            node_strengths_out[vertex_index(v,g)] += w
        end
    end
    return MGraph(g, node_strengths_in, node_strengths_out, zeros(num_vertices(g)), edge_weights, float(total_weight))
end

function collapse_graph{V}(mg::MGraph{V}, membership::Vector{Int})
    length(membership) == num_vertices(mg.graph) || error("wrong membership length")
    g = mg.graph
    nb_comm = maximum(membership)

    collapsed_node_strength_in = is_directed(g) ? zeros(nb_comm) : Float64[]
    collapsed_node_strength_out = zeros(nb_comm)
    collapsed_node_self_weights = zeros(nb_comm)
    collapsed_edge_weights = [Dict{Int,Float64}() for i=1:nb_comm]

    for e in edges(g)
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
        u = source(e,g)
        v = target(e,g)
        u_comm = membership[vertex_index(u,g)]
        v_comm = membership[vertex_index(v,g)]
<<<<<<< HEAD
        # for unidrected
        u_comm, v_comm = minmax(u_comm, v_comm)
        # excluding self-loop
        if u_comm != v_comm
            collapsed_strengths[u_comm] += w
            collapsed_strengths[v_comm] += w
            collapsed_edge_weights[u_comm][v_comm] = get(collapsed_edge_weights[u_comm], v_comm, 0.0) + w
        else
            collapsed_self_weights[u_comm] += 2w
        end
    end

    collapsed_graph = simple_graph(nb_comm, is_directed=false)
=======
        w = edge_property(mg.edge_weights, e, g)
        collapsed_node_strength_out[u_comm] += w
        if is_directed(g)
            collapsed_node_strength_in[v_comm] += w
        else
            collapsed_node_strength_out[v_comm] += w
            u_comm, v_comm = minmax(u_comm, v_comm)
        end

        if u_comm == v_comm
            # self loop weight
            collapsed_node_self_weights[u_comm] += w
        else
            collapsed_edge_weights[u_comm][v_comm] = get(collapsed_edge_weights[u_comm], v_comm, 0.0) + w
        end
    end

    collapsed_graph = simple_graph(nb_comm, is_directed=is_directed(g))
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
    collapsed_weights = Float64[]

    for u=1:nb_comm
        for (v,w) in collapsed_edge_weights[u]
            add_edge!(collapsed_graph, u, v)
            push!(collapsed_weights, w)
        end
    end

<<<<<<< HEAD
    MGraph(collapsed_graph, collapsed_strengths, collapsed_self_weights, VectorEdgePropertyInspector(collapsed_weights), mg.total_weight)
end

# type of flow partition
type MPartition{V}
    mgraph::MGraph{V}
    membership::Vector{Int}
    weights_inner::Vector{Float64}
    weights_out::Vector{Float64}
=======
    MGraph(collapsed_graph, collapsed_node_strength_in, collapsed_node_strength_out,
           collapsed_node_self_weights, VectorEdgePropertyInspector(collapsed_weights), mg.total_weight)
end

# type of modularity partition
type MPartition{V}
    mgraph::MGraph{V}
    membership::Vector{Int}
    total_weight_inner::Vector{Float64}
    total_weight_in::Vector{Float64}
    total_weight_out::Vector{Float64}
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
end

# construction
function mpartition{V}(g::AbstractGraph{V})
    n = num_vertices(g)
    mg = mgraph(g)
<<<<<<< HEAD
    MPartition(mg, collect(1:n), zeros(n), mg.strengths)
=======
    MPartition{V}(mg, collect(1:n), mg.node_self_weights, mg.node_strengths_in, mg.node_strengths_out)
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
end

function mpartition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T})
    n = num_vertices(g)
<<<<<<< HEAD
    mg = mgraph(g)
    MPartition(mg, collect(1:n), zeros(n), mg.strengths)
=======
    mg = mgraph(g, weights)
    MPartition{V}(mg, collect(1:n), mg.node_self_weights, mg.node_strengths_in, mg.node_strengths_out)
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
end

function collapse_partition{V}(mp::MPartition{V})
    collapsed_mg = collapse_graph(mp.mgraph, mp.membership)
    n = num_vertices(collapsed_mg.graph)
<<<<<<< HEAD
    MPartition(collapsed_mg, collect(1:n), collapsed_mg.self_weights, collapsed_mg.strengths)
end

function update_partition!{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    length(mp.membership) == num_vertices(g) || error("wrong membership length")


    nb_comm = maximum(mp.membership)
    mp.weights_out = zeros(nb_comm)
    mp.weights_inner = zeros(nb_comm)
    for e in edges(g)
        w = edge_property(mg.edge_weights, e, g)
        u = source(e,g)
        v = target(e,g)
        u_comm = mp.membership[vertex_index(u,g)]
        v_comm = mp.membership[vertex_index(v,g)]
        if u_comm != v_comm
            mp.weights_out[u_comm] += w
            mp.weights_out[v_comm] += w
        else
            mp.weights_inner[v_comm] += 2w
        end
    end
end



=======
    MPartition(collapsed_mg, mp.membership, collapsed_mg.node_self_weights,
               collapsed_mg.node_strengths_in, collapsed_mg.node_strengths_out)
end

>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
"Move a node to a new community and update the partition, this also removes any empty communities."
function move_node!{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph

    u_idx = vertex_index(u, g)
    old_comm = mp.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
<<<<<<< HEAD
        mp.weights_inner[old_comm] -= mg.self_weights[u_idx]
        mp.weights_inner[new_comm] += mg.self_weights[u_idx]
=======
        u_self_weight = mg.node_self_weights[u_idx]
        # change of inner weight of the old community for self loop
        mp.total_weight_inner[old_comm] -= u_self_weight

        # change of inner weight of the new community for self loop
        mp.total_weight_inner[new_comm] += u_self_weight

>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
        for e in out_edges(u, g)
            v = target(e, g)
            v_idx = vertex_index(v,g)
            v_comm = mp.membership[v_idx]
            w = edge_property(mg.edge_weights, e, g)
<<<<<<< HEAD
            # change of exit probability of the old community
            if v_comm != old_comm
                mp.weights_out[old_comm] -= w
            else
                mp.weights_out[old_comm] += w
                mp.weights_inner[old_comm] -= 2w
=======
            # change of out weight of the old community
            if v_comm != old_comm
                mp.total_weight_out[old_comm] -= w
            else
                mp.total_weight_inner[old_comm] -= w
                if is_directed(g)
                    mp.total_weight_in[old_comm] += w
                else
                    mp.total_weight_out[old_comm] += w
                end
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
            end
            # change of exit probability of the new community
            #if !in(v, partition.community[new_comm].nodes)
            if v_comm != new_comm
<<<<<<< HEAD
                mp.weights_out[new_comm] += w
            else
                mp.weights_out[new_comm] -= w
                mp.weights_inner[new_comm] += 2w
            end
        end

        # update the membership vector
        mp.membership[u_idx] = new_comm
    end
end

function renumber_communities!{V}(mp::MPartition{V})
    idx_map = Dict{Int,Int}()
    j = 1
    for i in mp.membership
        if !haskey(idx_map, i)
            idx_map[i] = j
            j += 1
        end
    end

    nb_comm = length(idx_map)
    new_weights_out = Array(Float64, nb_comm)
    new_weights_inner = Array(Float64, nb_comm)
    for (old_comm_idx, new_comm_idx) in idx_map
        new_weights_out[new_comm_idx] = mp.weights_out[old_comm_idx]
        new_weights_inner[new_comm_idx] = mp.weights_inner[old_comm_idx]
    end

    @inbounds for i=1:num_vertices(mp.mgraph.graph)
        mp.membership[i] = idx_map[mp.membership[i]]
    end
    mp.weights_out = new_weights_out
    mp.weights_inner = new_weights_inner
end

function from_coarser_partition!{V}(mp::MPartition{V}, coarser_mp::MPartition{V})
    for u=1:length(mp.membership)
        # what is the community of the node
        u_comm_level1 = mp.membership[u]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_mp.membership[u_comm_level1]
        mp.membership[u] = u_comm_level2
    end
    update_partition!(mp)
end

function diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u, g)
    old_comm = mp.membership[u_idx]
    diff_total = 0.0
    if new_comm != old_comm && length(unique(mp.membership)) > 2
        w_to_old = 0.0
        w_to_new = 0.0
        w_notto_old = 0.0
        w_notto_new = 0.0
        k_out = mg.strengths[u_idx]
        self_weight = mg.self_weights[u_idx]
        K_out_old = mp.weights_out[old_comm]
        K_out_new = mp.weights_out[new_comm]
        for e in out_edges(u,g)
            w = edge_property(mg.edge_weights, e, g)
            v = target(e, g)
            v_idx = vertex_index(v, g)
            v_comm = mp.membership[v_idx]
            if v_comm == old_comm
                w_to_old += w
            else
                w_notto_old += w
            end
            if v_comm == new_comm
                w_to_new += w
            else
                w_notto_new += w
            end
        end
        total_weight = 2*mg.total_weight
        diff_old = -(2w_to_old+self_weight) - ((w_to_old-w_notto_old)^2+2*(w_to_old-w_notto_old)*K_out_old)/total_weight
        diff_new = 2w_to_new + self_weight -((w_notto_new-w_to_new)^2+2*(w_notto_new-w_to_new)*K_out_new)/total_weight
        diff_total = (diff_new + diff_old)/total_weight
    end
    diff_total
end

function quality{V}(mp::MPartition{V})
    Q = 0.0
    total_weight = 2*mp.mgraph.total_weight
    if length(unique(mp.membership)) > 1
        for i=1:length(mp.weights_inner)
            Q += mp.weights_inner[i] - mp.weights_out[i]*mp.weights_out[i]/total_weight
        end
    end
    Q/total_weight
end

"get neighbor communities of node u"
function get_neigh_comms{V}(mp::MPartition{V}, u::V)
    g = mp.mgraph.graph
    neigh_comms = Set{Int}()
    for v in out_neighbors(u, g)
        v_idx = vertex_index(v, g)
        push!(neigh_comms, mp.membership[v_idx])
    end
    neigh_comms
end

"Move nodes to other communities depending on how other communities are
considered."
function move_nodes!{V}(mp::MPartition{V}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    g = mp.mgraph.graph
    memb = mp.membership
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
                max_diff = 0.0
                new_comm = u_comm


                # In which communities are its neighbours
                #neigh_comms = collect(get_neigh_comms(fp, u))
                neigh_comms = get_neigh_comms(mp, u)
                # tie break randomly, it has been found that this strategy can improv result
                #shuffle!(neigh_comms)
                # Loop through the communities of the neighbours
                for comm in neigh_comms
                    # Calculate the possible difference of the moving the node to that community
                    possible_diff = diff_move(mp, u, comm)
                    # We're only inserested in the maximum
                    if possible_diff > max_diff
                        max_diff = possible_diff
                        new_comm = comm
                    end
                end

                if new_comm != u_comm
                    # keep track of quality change
                    once_diff += max_diff
                    # actually move the node
                    move_node!(mp, u, new_comm)
                    num_moves += 1
                end
            end # if out_degree
        end # for

        # keep track of total difference over multiple loops
        total_diff += once_diff
    end # while
    renumber_communities!(mp)
    total_diff
end

"Optimize the provided partition."
function optimize_partition!{V}(mp::MPartition{V}; ϵ=1e-5, δ=1e-2, max_itr=10000)
    # do one iteration of optimisation
    once_diff = move_nodes!(mp, ϵ=ϵ, δ=δ, max_itr=max_itr)

    # as long as there remains imprvoment iterate
    while once_diff > ϵ
        # first create collapsed partition
        collapsed_mp = collapse_partition(mp)
        # Optimise partition for collapsed graph
        once_diff = move_nodes!(collapsed_mp, ϵ=ϵ, δ=δ, max_itr=max_itr)
        # Make sure improvement on coarser scale is reflected on the
        # scale of the graph as a whole.
        from_coarser_partition!(mp, collapsed_mp)
    end

    # We renumber the communities to make sure we stick in the range
    # 1,...,r for r communities.
    # By default, we number the communities in decreasing order of size,
    # so that 1 is the largest community, 2 the second largest, etc...
    renumber_communities!(mp)
    # Return the quality of the current partition.
    quality(mp)
end

g = simple_graph(4, is_directed=false)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 1)
add_edge!(g, 1, 4)
using GraphPlot
g = graphfamous("football")
mp = mpartition(g)
optimize_partition!(mp)
mp.membership
move_nodes!(mp)
from_coarser_partition!(mp, cmp)
mp.mgraph.graph
mp.membership
cmp = collapse_partition(mp)
cmp.weights_inner
mp.mgraph.graph
cmp.mgraph.total_weight
move_nodes!(cmp)
cmp.membership
cmp.mgraph.

mp.membership
diff_move(mp, 1, 4)
mp.weights_inner
mp.membership = [1,1,2,2]
update_partition!(mp)
move_node!(mp, 1, 4)
renumber_communities!(mp)
mp.weights_inner
quality(mp)
move_node!(mp, 1, 2)
quality(mp)
mp.weights_inner
diff_move(mp, 1, 2)
update_partition!(mp)
mp.weights_inner
cmp = collapse_partition(mp)
cmp.weights_out
mp = mpartition(g, [1.0,1.0,1.0,1.0])
mp.membership = [1,1,1,1]
cmp.weights_out
mp.weights_out
move_node!(mp, 1, 2)
mp.membership
renumber_communities!(mp)
mp.weights_out
mp.membership
from_coarser_partition!(mp, cmp)
mp.membership
=======
                mp.total_weight_out[new_comm] += w
            else
                mp.total_weight_inner[new_comm] += w
                if is_directed(g)
                    mp.total_weight_in[new_comm] -= w
                else
                    mp.total_weight_out[new_comm] -= w
                end
            end
        end

        # for directed graph only
        if is_directed(g)
            for e in in_edges(u, g)
                v = target(e, g)
                v_idx = vertex_index(v,g)
                v_comm = mp.membership[v_idx]
                w = edge_property(mg.edge_weights, e, g)
                # change of out weight of the old community
                if v_comm != old_comm
                    mp.total_weight_in[old_comm] -= w
                else
                    mp.total_weight_inner[old_comm] -= w
                    mp.total_weight_out[old_comm] += w
                end
                # change of exit probability of the new community
                #if !in(v, partition.community[new_comm].nodes)
                if v_comm != new_comm
                    mp.total_weight_in[new_comm] += w
                else
                    mp.total_weight_inner[new_comm] += w
                    mp.total_weight_out[new_comm] -= w
                end
            end
        end

        # update the membership vector
        mp.membership[u_idx] = new_comm
    end
end

mp = mpartition(g)
mp.membership
mp.total_weight_in
move_node!(mp, 2,1)
cmp = collapse_partition(mp)
cmp.mgraph.edge_weights
cg = collapse_graph(mp.mgraph, mp.membership)
g = simple_graph(4, is_directed=true)
add_edge!(g,1,2)
add_edge!(g,2,3)
add_edge!(g,3,1)
add_edge!(g,1,4)
mg = mgraph(g)
cmg = collapse_graph(mg, [1,1,2,2])
cmg.node_self_weights
>>>>>>> 243a36b67120c0a1ec41fd4a932bb07bf5556634
