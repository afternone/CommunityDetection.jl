"Move nodes to other communities depending on how other communities are
considered."
function move_nodes!{V}(partition::AbstractPartition{V}, ϵ=1e-5; δ=1e-2,
                       max_itr=10000, random_order=true, consider_comms=:all_neigh_comms)
    g = graph(partition)
    memb = membership(partition)
    communities = community(partition)
    itr = 0 # number of iterations
    total_diff = 0.0 # total difference while moving nodes
    once_diff = 2ϵ # difference for one loop
    n = num_vertices(g)
    num_moves = 2n # number of moved nodes during one loop
    # As long as we keep on improving and we don't exceed the
    # maximum number of iterations and number of moves.
    while once_diff > ϵ && num_moves > n*δ && itr < max_itr
        # increase number of iterations
        itr += 1

        # initialize number of moves and imporvment
        num_moves = 0
        once_diff = 0.0

        # establish vertex order
        # we normally initialize the normal vertex order
        # of considering node 0,1,...
        vertex_order = collect(vertices(g))

        # But if we use a random order, we shuffle this order.
        if random_order
            shuffle!(vertex_order)
        end

        # for each node
        for u in vertex_order
            # only take into account nodes of degree higher than zero
            if out_degree(u, g) > 0
                u_idx = vertex_index(u, g)
                # What is the current community of the node
                u_comm = memb[u_idx]
                # What is the difference per community if we move the node to one of
                # the other communities, and what is the extreme difference?
                extreme_diff = 0.0
                extreme_comm = u_comm

                # Keep track of the possible improvements and (neighbouring) communities.
                if consider_comms == :all_comms
                    # loop through all communities
                    for comm in communities
                        # Calculate the possible improvement of the moving the node to that community
                        possible_diff = diff_move(partition, u, comm)
                        # We're only inserested in the extreme
                        # TO DO, tie break randomly
                        if is_right_direction(partition, extreme_diff, possible_diff)
                            extreme_diff = possible_diff
                            extreme_comm = comm
                        end
                    end
                elseif consider_comms == :all_neigh_comms
                    # In which communities are its neighbours
                    #neigh_comms = get_neigh_comms(partition, u)
                    # tie break randomly, it has been found that this strategy can improv result
                    neigh_comms = collect(get_neigh_comms(partition, u))
                    shuffle!(neigh_comms)
                    # Loop through the communities of the neighbours
                    for comm in neigh_comms
                        # Calculate the possible difference of the moving the node to that community
                        possible_diff = diff_move(partition, u, comm)
                        # We're only inserested in the extreme
                        if is_right_direction(partition, extreme_diff, possible_diff)
                            extreme_diff = possible_diff
                            extreme_comm = comm
                        end
                    end
                elseif consider_comms == :rand_comm
                    # Select a random community
                    neigh_comm = memb[rand(1:n)]
                    possible_diff = diff_move(partition, u, neigh_comm)
                    if is_right_direction(partition, 0.0, possible_diff)
                        extreme_diff = possible_diff
                        extreme_comm = neigh_comm
                    end
                elseif consider_comms == :rand_neigh_comm
                    u_neighbors = out_neighbors(u, g)
                    u_neighbors_idx = Int[vertex_index(v, g) for v in u_neighbors]
                    neigh_comm = memb[u_neighbors_idx[rand(1:length(u_neighbors_idx))]]
                    possible_diff = diff_move(partition, u, neigh_comm)
                    if is_right_direction(partition, 0.0, possible_diff)
                        extreme_diff = possible_diff
                        extreme_comm = neigh_comm
                    end
                end # if consider_comms

                if extreme_comm != u_comm
                    # keep track of quality change
                    once_diff += extreme_diff
                    # actually move the node
                    move_node!(partition, u, extreme_comm)
                    num_moves += 1
                end

            end # if out_degree
        end # for

        # keep track of total difference over multiple loops
        total_diff += once_diff
    end # while
    renumber_communities!(partition)
    total_diff
end

"Optimize the provided partition."
function optimize_partition!{V}(partition::AbstractPartition{V}, ϵ=1e-5; kargs...)
    # do one iteration of optimisation
    once_diff = move_nodes!(partition, ϵ; kargs...)
    # as long as there remains imprvoment iterate
    while still_running(partition, once_diff, ϵ)
        # first create collapsed partition
        collapsed_partition = collapse_partition(partition)
        # Optimise partition for collapsed graph
        once_diff = move_nodes!(collapsed_partition, ϵ; kargs...)
        # Make sure improvement on coarser scale is reflected on the
        # scale of the graph as a whole.
        from_coarser_partition!(partition, collapsed_partition)
    end

    # We renumber the communities to make sure we stick in the range
    # 1,...,r for r communities.
    # By default, we number the communities in decreasing order of size,
    # so that 1 is the largest community, 2 the second largest, etc...
    renumber_communities!(partition)
    # Return the quality of the current partition.
    quality(partition)
end

function find_partition!{V}(partition::AbstractPartition{V}, ϵ=1e-5; kargs...)
    quality1 = optimize_partition!(partition, ϵ; kargs...)
    quality2 = optimize_partition!(partition, ϵ; kargs...)

    while abs(quality2-quality1) > ϵ
        quality1 = quality2
        quality2 = optimize_partition!(partition, ϵ; kargs...)
    end
    quality2
end

"Move nodes to other communities depending on how other communities are
considered."
function move_nodes!{V}(partition::MPartition{V}, _diff_move::Function = modularity_diff_move, ϵ=1e-5; δ=1e-2,
                       max_itr=10000, random_order=true, consider_comms=:all_neigh_comms)
    g = graph(partition)
    memb = membership(partition)
    communities = community(partition)
    itr = 0 # number of iterations
    total_diff = 0.0 # total difference while moving nodes
    once_diff = 2ϵ # difference for one loop
    n = num_vertices(g)
    num_moves = 2n # number of moved nodes during one loop
    # As long as we keep on improving and we don't exceed the
    # maximum number of iterations and number of moves.
    while once_diff > ϵ && num_moves > n*δ && itr < max_itr
        # increase number of iterations
        itr += 1

        # initialize number of moves and imporvment
        num_moves = 0
        once_diff = 0.0

        # establish vertex order
        # we normally initialize the normal vertex order
        # of considering node 0,1,...
        vertex_order = collect(vertices(g))

        # But if we use a random order, we shuffle this order.
        if random_order
            shuffle!(vertex_order)
        end

        # for each node
        for u in vertex_order
            # only take into account nodes of degree higher than zero
            if out_degree(u, g) > 0
                u_idx = vertex_index(u, g)
                # What is the current community of the node
                u_comm = memb[u_idx]
                # What is the difference per community if we move the node to one of
                # the other communities, and what is the extreme difference?
                extreme_diff = 0.0
                extreme_comm = u_comm

                # Keep track of the possible improvements and (neighbouring) communities.
                if consider_comms == :all_comms
                    # loop through all communities
                    for comm in communities
                        # Calculate the possible improvement of the moving the node to that community
                        possible_diff = _diff_move(partition, u, comm)
                        # We're only inserested in the extreme
                        # TO DO, tie break randomly
                        if is_right_direction(partition, extreme_diff, possible_diff)
                            extreme_diff = possible_diff
                            extreme_comm = comm
                        end
                    end
                elseif consider_comms == :all_neigh_comms
                    # In which communities are its neighbours
                    #neigh_comms = get_neigh_comms(partition, u)
                    # tie break randomly, it has been found that this strategy can improv result
                    neigh_comms = collect(get_neigh_comms(partition, u))
                    shuffle!(neigh_comms)
                    # Loop through the communities of the neighbours
                    for comm in neigh_comms
                        # Calculate the possible difference of the moving the node to that community
                        possible_diff = _diff_move(partition, u, comm)
                        # We're only inserested in the extreme
                        if is_right_direction(partition, extreme_diff, possible_diff)
                            extreme_diff = possible_diff
                            extreme_comm = comm
                        end
                    end
                elseif consider_comms == :rand_comm
                    # Select a random community
                    neigh_comm = memb[rand(1:n)]
                    possible_diff = _diff_move(partition, u, neigh_comm)
                    if is_right_direction(partition, 0.0, possible_diff)
                        extreme_diff = possible_diff
                        extreme_comm = neigh_comm
                    end
                elseif consider_comms == :rand_neigh_comm
                    u_neighbors = out_neighbors(u, g)
                    u_neighbors_idx = Int[vertex_index(v, g) for v in u_neighbors]
                    neigh_comm = memb[u_neighbors_idx[rand(1:length(u_neighbors_idx))]]
                    possible_diff = _diff_move(partition, u, neigh_comm)
                    if is_right_direction(partition, 0.0, possible_diff)
                        extreme_diff = possible_diff
                        extreme_comm = neigh_comm
                    end
                end # if consider_comms

                if extreme_comm != u_comm
                    # keep track of quality change
                    once_diff += extreme_diff
                    # actually move the node
                    move_node!(partition, u, extreme_comm)
                    num_moves += 1
                end

            end # if out_degree
        end # for

        # keep track of total difference over multiple loops
        total_diff += once_diff
    end # while
    renumber_communities!(partition)
    total_diff
end

"Optimize the provided partition."
function optimize_partition!{V}(partition::MPartition{V}, ϵ=1e-5; δ=1e-2, method = :Modularity,
                       max_itr=10000, random_order=true, consider_comms=:all_neigh_comms, resolution_parameter::Float64=0.1)
    if method == :Modularity
        _diff_move = modularity_diff_move
        _quality = modularity_quality
    elseif method == :CPM
        _diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int) = CPM_diff_move(mp, u, new_comm, resolution_parameter)
        _quality{V}(mp::MPartition{V}) = CPM_quality(mp, resolution_parameter)
    elseif method == :Significance
        if !all(x->x==one(eltype(partition.mgraph.edge_weights)), partition.mgraph.edge_weights)
            error("Significance is not suited for optimisation on weighted graphs. Please consider a different method.")
        end
        _diff_move = significance_diff_move
        _quality = significance_quality
    elseif method == :Surprise
        _diff_move = surprise_diff_move
        _quality = surprise_quality
    elseif method == :RBConfiguration
        _diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int) = RBC_diff_move(mp, u, new_comm, resolution_parameter)
        _quality{V}(mp::MPartition{V}) = RBC_quality(mp, resolution_parameter)
    elseif method == :RBER
        _diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int) = RBER_diff_move(mp, u, new_comm, resolution_parameter)
        _quality{V}(mp::MPartition{V}) = RBER_quality(mp, resolution_parameter)
    else
        error("Non-existing method for optimization specified.")
    end

    # do one iteration of optimisation
    once_diff = move_nodes!(partition, _diff_move, ϵ, δ=δ, max_itr=max_itr, random_order=random_order, consider_comms=consider_comms)
    # as long as there remains imprvoment iterate
    while once_diff > ϵ
        # first create collapsed partition
        collapsed_partition = collapse_partition(partition)
        # Optimise partition for collapsed graph
        once_diff = move_nodes!(collapsed_partition, _diff_move, ϵ, δ=δ, max_itr=max_itr, random_order=random_order, consider_comms=consider_comms)
        # Make sure improvement on coarser scale is reflected on the
        # scale of the graph as a whole.
        from_coarser_partition!(partition, collapsed_partition)
    end

    # We renumber the communities to make sure we stick in the range
    # 1,...,r for r communities.
    # By default, we number the communities in decreasing order of size,
    # so that 1 is the largest community, 2 the second largest, etc...
    renumber_communities!(partition)
    # Return the quality of the current partition.
    _quality(partition)
end
