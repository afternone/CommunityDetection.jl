# type defined
type MGroup{V}
    nodes::Set{V}
    csize::Int
    weight_inner::Float64
    weight_in::Float64
    weight_out::Float64
end

type MPartition{V} <: AbstractPartition{V}
    mgraph::MGraph{V}
    membership::Vector{Int}
    community::Dict{Int,MGroup{V}}
    total_weight_in_all_comms::Float64
    total_possible_edges_in_all_comms::Int
end

# require interface
graph{V}(mp::MPartition{V}) = mp.mgraph.graph
membership{V}(mp::MPartition{V}) = mp.membership
membership{V}(mp::MPartition{V}, u::V) = mp.membership[vertex_index(u, mp.mgraph.graph)]
community{V}(mp::MPartition{V}) = keys(mp.community)
is_right_direction{V}(mp::MPartition{V}, quality, new_quality) = new_quality > quality
still_running{V}(mp::MPartition{V}, once_diff, ϵ) = once_diff > ϵ

# construction
function mpartition{V,T<:Real}(g::AbstractGraph{V}, edge_weights::Vector{T}=ones(num_edges(g)),
                               node_sizes::Vector{Int}=ones(Int,num_vertices(g)), correct_self_loops::Bool=false)
    mg = mgraph(g, edge_weights, node_sizes, correct_self_loops)
    mpartition(mg)
end

function mpartition{V}(mg::MGraph{V})
    n = num_vertices(mg.graph)
    mp = MPartition{V}(mg, collect(1:n), Dict{Int,MGroup{V}}(), 0.0, 0)
    update_partition!(mp)
    mp
end

# mutation
function update_partition!{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    n = num_vertices(g)
    nb_comms = maximum(mp.membership)
    # reset partition
    node_sets = Array(Set{V}, nb_comms)
    for i=1:nb_comms
        node_sets[i] = Set{V}()
    end

    w_inner = zeros(nb_comms)
    w_in = zeros(nb_comms)
    w_out = zeros(nb_comms)
    csize = zeros(Int, nb_comms)

    mp.total_weight_in_all_comms = 0.0
    for u in vertices(g)
        u_idx = vertex_index(u,g)
        u_comm = mp.membership[u_idx]
        # add this node to the communtiy sets
        push!(node_sets[u_comm], u)
        # update the communtiy size
        csize[u_comm] += mg.node_sizes[u_idx]

        # loop over all out edges
        for e in out_edges(u,g)
            e_idx = edge_index(e,g)
            v = target(e,g)
            v_idx = vertex_index(v,g)
            v_comm = mp.membership[v_idx]
            # get the weight of the edge
            w = mg.edge_weights[e_idx]
            # Add weight to the outgoing weight of community of u
            w_out[u_comm] += w
            # Add weight to the incoming weight of community of v
            w_in[v_comm] += w
            # if it is an edge within a community
            if u_comm == v_comm
                if !is_directed(g)
                    w /= 2
                end
                w_inner[u_comm] += w
                mp.total_weight_in_all_comms += w
            end
        end
    end
    mp.total_possible_edges_in_all_comms = 0
    for c=1:nb_comms
        n_c = csize[c]
        possible_edges = 0
        if mg.correct_self_loops
            possible_edges = round(Int, n_c*n_c/(2.0 - Int(is_directed(g))))
        else
            possible_edges = round(Int, n_c*(n_c-1)/(2.0 - Int(is_directed(g))))
        end
        mp.total_possible_edges_in_all_comms += possible_edges
    end

    empty!(mp.community)
    for i=1:nb_comms
        if !isempty(node_sets[i])
            mp.community[i] = MGroup{V}(node_sets[i], csize[i], w_inner[i], w_in[i], w_out[i])
        end
    end
end

function move_node!{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    node_size = mg.node_sizes[u_idx]
    old_comm = mp.membership[u_idx]
    old_csize = mp.community[old_comm].csize
    new_csize = mp.community[new_comm].csize
    mp.total_possible_edges_in_all_comms += 2.0*node_size*(new_csize - old_csize + node_size)/(2.0 - Float64(is_directed(g)))

    # remove from old community
    delete!(mp.community[old_comm].nodes, u)
    mp.community[old_comm].csize -= node_size

    # add to new community
    push!(mp.community[new_comm].nodes, u)
    mp.community[new_comm].csize += node_size

    for e in out_edges(u, g)
        e_idx = edge_index(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_comm = mp.membership[v_idx]
        w = mg.edge_weights[e_idx]
        mp.community[old_comm].weight_out -= w
        mp.community[new_comm].weight_out += w
        int_weight = w/(is_directed(g) ? 1.0 : 2.0)/( u == v ? 2.0 : 1.0)
        if old_comm == v_comm
            # remove the internal weight
            mp.community[old_comm].weight_inner -= int_weight
            mp.total_weight_in_all_comms -= int_weight
        end
        if new_comm == v_comm || u == v
            # add the internal weight
            mp.community[new_comm].weight_inner += int_weight
            mp.total_weight_in_all_comms += int_weight
        end
    end

    for e in in_edges(u, g)
        e_idx = edge_index(e,g)
        v = source(e,g)
        v_idx = vertex_index(v,g)
        v_comm = mp.membership[v_idx]
        w = mg.edge_weights[e_idx]
        mp.community[old_comm].weight_in -= w
        mp.community[new_comm].weight_in += w
        int_weight = w/(is_directed(g) ? 1.0 : 2.0)/( u == v ? 2.0 : 1.0)
        if old_comm == v_comm
            # remove the internal weight
            mp.community[old_comm].weight_inner -= int_weight
            mp.total_weight_in_all_comms -= int_weight
        end
        if new_comm == v_comm || u == v
            # add the internal weight
            mp.community[new_comm].weight_inner += int_weight
            mp.total_weight_in_all_comms += int_weight
        end
    end

    # if the old community is empty after remove node u, we remove it
    if isempty(mp.community[old_comm].nodes)
        delete!(mp.community, old_comm)
    end

    # update the membership vector
    mp.membership[u_idx] = new_comm
end

function collapse_partition{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    n = num_vertices(g)
    m = num_edges(g)
    nb_comm = length(mp.community)

    collapsed_edge_weights = Array(Dict{Int,Float64}, nb_comm)
    for i=1:nb_comm
        collapsed_edge_weights[i] = Dict{Int,Float64}()
    end

    for e in edges(g)
        e_idx = edge_index(e,g)
        w = mg.edge_weights[e_idx]
        u = source(e,g)
        v = target(e,g)
        u_idx = vertex_index(u,g)
        v_idx = vertex_index(v,g)
        u_comm = mp.membership[u_idx]
        v_comm = mp.membership[v_idx]

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

    # Now create vector for edges, first determined the number of edges
    n_collapsed = nb_comm
    m_collapsed = 0

    for itr in collapsed_edge_weights
        m_collapsed += length(itr)
    end

    collapsed_graph = simple_graph(n_collapsed, is_directed=is_directed(g))

    collapsed_weights = Float64[]
    total_collapsed_weight = 0.0

    for u=1:n_collapsed
        for (v,w) in collapsed_edge_weights[u]
            add_edge!(collapsed_graph, u, v)
            push!(collapsed_weights, w)
            total_collapsed_weight += w
        end
    end

    if abs(total_collapsed_weight-mg.total_weight) > 1.0e-6
        error("Total collapsed weight is not equal to original weight.")
    end

    if num_vertices(collapsed_graph) != nb_comm
        error("Something went wrong with collapsing the graph.")
    end

    # calculate new node sizes
    csizes = zeros(Int, n_collapsed)
    for (i,c) in mp.community
        csizes[i] = c.csize
    end

    collapsed_mg = mgraph(collapsed_graph, collapsed_weights, csizes, mg.correct_self_loops)
    mpartition(collapsed_mg)
end

function diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_community(mg, u, old_comm, mp.membership)
        w_from_old = weight_from_community(mg, u, old_comm, mp.membership)
        w_to_new = weight_to_community(mg, u, new_comm, mp.membership)
        w_from_new = weight_from_community(mg, u, new_comm, mp.membership)
        k_out = mg.strength_out[u_idx]
        k_in = mg.strength_in[u_idx]
        self_weight = mg.node_self_weights[u_idx]
        K_out_old = mp.community[old_comm].weight_out
        K_in_old = mp.community[old_comm].weight_in
        K_out_new = mp.community[new_comm].weight_out + k_out
        K_in_new = mp.community[new_comm].weight_in + k_in
        total_weight = mg.total_weight*(2.0 - Float64(is_directed(g)))
        diff_old = (w_to_old - k_out*K_in_old/total_weight) + (w_from_old - k_in*K_out_old/total_weight)
        diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) +
            (w_from_new + self_weight - k_in*K_out_new/total_weight)
        diff = (diff_new - diff_old)/total_weight
    end
    diff
end

function modularity_diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_community(mg, u, old_comm, mp.membership)
        w_from_old = weight_from_community(mg, u, old_comm, mp.membership)
        w_to_new = weight_to_community(mg, u, new_comm, mp.membership)
        w_from_new = weight_from_community(mg, u, new_comm, mp.membership)
        k_out = mg.strength_out[u_idx]
        k_in = mg.strength_in[u_idx]
        self_weight = mg.node_self_weights[u_idx]
        K_out_old = mp.community[old_comm].weight_out
        K_in_old = mp.community[old_comm].weight_in
        K_out_new = mp.community[new_comm].weight_out + k_out
        K_in_new = mp.community[new_comm].weight_in + k_in
        total_weight = mg.total_weight*(2.0 - Float64(is_directed(g)))
        diff_old = (w_to_old - k_out*K_in_old/total_weight) + (w_from_old - k_in*K_out_old/total_weight)
        diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) +
            (w_from_new + self_weight - k_in*K_out_new/total_weight)
        diff = (diff_new - diff_old)/total_weight
    end
    diff
end

function quality{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    total_weight = mg.total_weight*(2.0 - Float64(is_directed(g)))
    Q = 0.0
    for c in values(mp.community)
        w = c.weight_inner
        w_out = c.weight_out
        w_in = c.weight_in
        Q += w - w_out*w_in/((is_directed(g) ? 1.0 : 4.0)*mg.total_weight)
    end
    (2.0 - Float64(is_directed(g)))*Q / total_weight
end

function modularity_quality{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    total_weight = mg.total_weight*(2.0 - Float64(is_directed(g)))
    Q = 0.0
    for c in values(mp.community)
        w = c.weight_inner
        w_out = c.weight_out
        w_in = c.weight_in
        Q += w - w_out*w_in/((is_directed(g) ? 1.0 : 4.0)*mg.total_weight)
    end
    (2.0 - Float64(is_directed(g)))*Q / total_weight
end

# surprise partition
function surprise_diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    nsize = mp.community[old_comm].csize
    diff = 0.0
    if new_comm != old_comm
        normalise = 2.0 - Float64(is_directed(g))
        m = mg.total_weight
        n = mg.total_size
        n2 = mg.correct_self_loops ? n*n/normalise : n*(n-1)/normalise
        # before move
        mc = mp.total_weight_in_all_comms
        nc2 = mp.total_possible_edges_in_all_comms

        # to old comm
        n_old = mp.community[old_comm].csize
        sw = mg.node_self_weights[u_idx]
        wtc = weight_to_community(mg, u, old_comm, mp.membership) - sw
        wfc = weight_from_community(mg, u, old_comm, mp.membership) - sw
        m_old = wtc/normalise + wfc/normalise + sw

        # to new comm
        n_new = mp.community[new_comm].csize
        wtc = weight_to_community(mg, u, new_comm, mp.membership)
        wfc = weight_from_community(mg, u, new_comm, mp.membership)
        sw = mg.node_self_weights[u_idx]
        m_new = wtc/normalise + wfc/normalise + sw

        q = mc/m
        s = nc2/n2
        q_new = (mc - m_old + m_new)/m
        s_new = (nc2 + 2nsize*(n_new - n_old + nsize)/normalise)/n2
        diff = m*(KL(q_new, s_new) - KL(q, s))
    end
    diff
end

function surprise_quality{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    normalise = 2.0 - Float64(is_directed(g))
    mc = mp.total_weight_in_all_comms
    nc2 = mp.total_possible_edges_in_all_comms
    m = mg.total_weight
    n = mg.total_size
    n2 = mg.correct_self_loops ? n*n/normalise : n*(n-1)/normalise
    q = mc/m
    s = nc2/n2
    S = m*KL(q,s)
    S
end

# CPM partition
function CPM_diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int, resolution_parameter::Float64)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_community(mg, u, old_comm, mp.membership)
        w_to_new = weight_to_community(mg, u, new_comm, mp.membership)
        w_from_old = weight_from_community(mg, u, old_comm, mp.membership)
        w_from_new = weight_from_community(mg, u, new_comm, mp.membership)
        nsize = mg.node_sizes[u_idx]
        csize_old = mp.community[old_comm].csize
        csize_new = mp.community[new_comm].csize
        self_weight = mg.node_self_weights[u_idx]
        possible_edge_difference_old = 0.0
        if mg.correct_self_loops
            possible_edge_difference_old = nsize*(2.0*csize_old - nsize)
        else
            possible_edge_difference_old = nsize*(2.0*csize_old - nsize - 1.0)
        end
        diff_old = w_to_old + w_from_old -
            self_weight - resolution_parameter*possible_edge_difference_old
        possible_edge_difference_new = 0.0
        if mg.correct_self_loops
             possible_edge_difference_new = nsize*(2.0*csize_new + nsize)
        else
             possible_edge_difference_new = nsize*(2.0*csize_new + nsize - 1.0)
        end
        diff_new = w_to_new + w_from_new + self_weight -
            resolution_parameter*possible_edge_difference_new
        diff = diff_new - diff_old
    end
    diff
end

function CPM_quality{V}(mp::MPartition{V}, resolution_parameter::Float64)
    mg = mp.mgraph
    g = mg.graph
    Q = 0.0
    for grp in values(mp.community)
        csize = grp.csize
        w = grp.weight_inner
        comm_possible_edges = mg.correct_self_loops ? csize*csize : csize*(csize-1)
        if !is_directed(g)
            comm_possible_edges /= 2
        end
        Q += w - resolution_parameter*comm_possible_edges
    end
    (2.0 - Float64(is_directed(g)))*Q
end

# RBConfiguration partition
function RBC_diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int, resolution_parameter::Float64)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_community(mg, u, old_comm, mp.membership)
        w_from_old = weight_from_community(mg, u, old_comm, mp.membership)
        w_to_new = weight_to_community(mg, u, new_comm, mp.membership)
        w_from_new = weight_from_community(mg, u, new_comm, mp.membership)
        k_out = mg.strength_out[u_idx]
        k_in = mg.strength_in[u_idx]
        self_weight = mg.node_self_weights[u_idx]
        K_out_old = mp.community[old_comm].weight_out
        K_in_old = mp.community[old_comm].weight_in
        K_out_new = mp.community[new_comm].weight_out
        K_in_new = mp.community[new_comm].weight_in
        total_weight = mg.total_weight*(2-Float64(is_directed(g)))
        diff_old = (w_to_old - resolution_parameter*k_out*K_in_old/total_weight) +
            (w_from_old - resolution_parameter*k_in*K_out_old/total_weight)
        diff_new = (w_to_new + self_weight - resolution_parameter*k_out*K_in_new/total_weight) +
            (w_from_new + self_weight - resolution_parameter*k_in*K_out_new/total_weight)
        diff = diff_new - diff_old
    end
    diff
end

function RBC_quality{V}(mp::MPartition{V}, resolution_parameter::Float64)
    mg = mp.mgraph
    g = mg.graph
    Q = 0.0
    for c in values(mp.community)
        w = c.weight_inner
        w_out = c.weight_out
        w_in = c.weight_in
        Q += w - resolution_parameter*w_out*w_in/((is_directed(g) ? 1.0 : 4.0)*mg.total_weight)
    end
    (2.0 - Float64(is_directed(g)))*Q
end

# RBER partition
function RBER_diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int, resolution_parameter::Float64)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_community(mg, u, old_comm, mp.membership)
        w_to_new = weight_to_community(mg, u, new_comm, mp.membership)
        w_from_old = weight_from_community(mg, u, old_comm, mp.membership)
        w_from_new = weight_from_community(mg, u, new_comm, mp.membership)
        nsize = mg.node_sizes[u_idx]
        csize_old = mp.community[old_comm].csize
        csize_new = mp.community[new_comm].csize
        self_weight = mg.node_self_weights[u_idx]
        possible_edge_difference_old = 0.0
        if mg.correct_self_loops
            possible_edge_difference_old = nsize*(2.0*csize_old - nsize)
        else
            possible_edge_difference_old = nsize*(2.0*csize_old - nsize - 1.0)
        end
        diff_old = w_to_old + w_from_old - self_weight -
            resolution_parameter*mg.density*possible_edge_difference_old
        possible_edge_difference_new = 0.0
        if mg.correct_self_loops
             possible_edge_difference_new = nsize*(2.0*csize_new + nsize)
        else
             possible_edge_difference_new = nsize*(2.0*csize_new + nsize - 1.0)
        end
        diff_new = w_to_new + w_from_new + self_weight -
            resolution_parameter*mg.density*possible_edge_difference_new
        diff = diff_new - diff_old
    end
    diff
end

function RBER_quality{V}(mp::MPartition{V}, resolution_parameter::Float64)
    mg = mp.mgraph
    g = mg.graph
    Q = 0.0
    for c in values(mp.community)
        csize = c.csize
        w = c.weight_inner
        comm_possible_edges = mg.correct_self_loops ? csize*csize : csize*(csize-1)
        if !is_directed(g)
            comm_possible_edges /= 2
        end
        Q += w - resolution_parameter*mg.density*comm_possible_edges
    end
    (2.0 - Float64(is_directed(g)))*Q
end

# Significance partition
function significance_diff_move{V}(mp::MPartition{V}, u::V, new_comm)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    nsize = mg.node_sizes[u_idx]
    diff = 0.0
    if new_comm != old_comm
        normalise = 2.0 - Float64(is_directed(g))
        p = mg.density
        # old community
        n_old = mp.community[old_comm].csize
        m_old = mp.community[old_comm].weight_inner
        q_old = 0.0
        if n_old > 1
            q_old = m_old/(n_old*(n_old-1)/normalise)
        end
        # old community after move
        n_oldx = n_old - nsize
        sw = mg.node_self_weights[u_idx]
        wtc = weight_to_community(mg, u, old_comm, mp.membership) - sw
        wfc = weight_from_community(mg, u, old_comm, mp.membership) - sw
        m_oldx = m_old - wtc/normalise - wfc/normalise - sw
        q_oldx = 0.0
        if n_oldx > 1
            q_oldx = m_oldx/(n_oldx*(n_oldx - 1)/normalise)
        end
        # new community
        n_new = mp.community[new_comm].csize
        m_new = mp.community[new_comm].weight_inner
        q_new = 0.0
        if n_new > 1
            q_new = m_new/(n_new*(n_new - 1)/normalise)
        end
        # new community after move
        n_newx = n_new + nsize
        wtc = weight_to_community(mg, u, new_comm, mp.membership)
        wfc = weight_from_community(mg, u, new_comm, mp.membership)
        sw = mg.node_self_weights[u_idx]
        m_newx = m_new + wtc/normalise + wfc/normalise + sw
        q_newx = 0.0
        if n_newx > 1
            q_newx = m_newx/(n_newx*(n_newx - 1)/normalise)
        end
        # calculate actual diff
        diff = - n_old*(n_old-1)*KL(q_old, p) +
            n_oldx*(n_oldx-1)*KL(q_oldx, p) -
            n_new*(n_new-1)*KL(q_new, p)+
            n_newx*(n_newx-1)*KL(q_newx, p)
    end
    diff
end

function significance_quality{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    S = 0.0
    p = mg.density
    for c in values(mp.community)
        n_c = c.csize
        m_c = c.weight_inner
        p_c = 0.0
        if n_c > 1
            p_c = m_c/(n_c*(n_c - 1)/(2.0 - Float64(is_directed(g))))
            S += KL(p_c, p)*n_c*(n_c - 1)
        end
    end
    S
end
