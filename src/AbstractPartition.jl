# the root type of all partitions
abstract AbstractPartition{V}

"Renumber the communities so that they are numbered 0,...,q-1 where q is the number of communities."
function renumber_communities!{V}(partition::AbstractPartition{V})
    g = graph(partition)
    csizes = Int[length(partition.community[i].nodes) for i in keys(partition.community)]
    perm_idx = sortperm(csizes, rev=true)
    temp_comms = collect(values(partition.community))
    empty!(partition.community)
    for (i, idx) in enumerate(perm_idx)
        partition.community[i] = temp_comms[idx]
    end

    for (i, group) in partition.community
        for u in group.nodes
            u_idx = vertex_index(u, g)
            partition.membership[u_idx] = i
        end
    end
end

function from_coarser_partition!{V}(partition::AbstractPartition{V}, coarser_partition::AbstractPartition{V})
    g = graph(partition)
    for u in vertices(g)
        u_idx = vertex_index(u, g)
        # what is the community of the node
        u_comm_level1 = partition.membership[u_idx]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_partition.membership[u_comm_level1]
        partition.membership[u_idx] = u_comm_level2
    end
    update_partition!(partition)
end

"Read new partition from another partition."
function from_partition!{V}(partition::AbstractPartition{V}, other_partition::AbstractPartition{V})
    #Assign the membership of every node in the supplied partition to the one in this partition
    partition.membership[:] = other_partition.membership[:]
    update_partition!(partition)
end

"get neighbor communities of node u"
function get_neigh_comms{V}(partition::AbstractPartition{V}, u::V)
    g = graph(partition)
    neigh_comms = Set{Int}()
    for v in out_neighbors(u, g)
        v_idx = vertex_index(v, g)
        push!(neigh_comms, partition.membership[v_idx])
    end
    neigh_comms
end
