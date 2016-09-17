using StatsBase

"""create GN benchmark graph with given expected value of intra-partition edges"""
function GN_network(k_in)
    #calculate edge probabilities
    k_out = 16-k_in
    p = k_in/31.0
    q = k_out/96
    #create graph/partitions with nodes
    g = simple_graph(128, is_directed=false)
    for k=1:4,i=1:31,j=i+1:32
        if rand() < p
            add_edge!(g, 32*(k-1)+i, 32*(k-1)+j)
        end
    end
    for k1=1:4,k2=k1+1:4
        if (k1!=k2)
            for i=1:32,j=1:32
                if rand()<q
                    add_edge!(g, 32*(k1-1)+i,32*(k2-1)+j)
                end
            end
        end
    end
    g
end

"""
Evaluate a partition of GN synthetic test graph, return fvcc score.

From Newman's article: (Fast algorithm for detecting community structure in networks)

    The criterion for deciding correct classification is as follows.
    We find the largest set of vertices that are grouped together by
    the algorithm in each of the four known communities. If the algorithm
    puts two or more of these sets in the same group, then all vertices
    in those sets are considered incorrectly classi- fied. Otherwise,
    they are considered correctly classified. All other vertices not in
    the largest sets are considered incorrectly classified.
"""
function fvcc(partition)
    # generate true partition
    real_part = [fill(1,32);fill(2,32);fill(3,32);fill(4,32)]
    fvcc(partition, real_part)
end

function fvcc(partition, real_partition)
    # get groups in the original clusters
    real_clusters = Dict{Int,Vector{Int}}()
    for (cluster,i) in zip(real_partition, 1:length(real_partition))
        if haskey(real_clusters, cluster)
            push!(real_clusters[cluster], i)
        else
            real_clusters[cluster] = [i]
        end
    end

    # get largest groups for each original group
    groups = Dict{Int, Pair{Int,Int}}()
    for (cluster,nodes) in real_clusters
        groups[cluster] = get_largest_group([partition[i] for i in nodes])
    end

    # add up sizes and check overlapping
    correct_count = 0
    for i in keys(real_clusters)
        # check for overlapping
        correct = true
        for j in keys(real_clusters)
            if i!=j && groups[i][1]==groups[j][1]
                correct = false
            end
        end
        # sum if correct
        if correct
            correct_count += groups[i][2]
        end
    end
    correct_count/length(partition)
end

""" Calculate largest group and its count in a list."""
function get_largest_group(input_list)
    group_counts = countmap(input_list)
    maxkey, maxvalue = next(group_counts, start(group_counts))[1]
    for (key, value) in group_counts
        if value > maxvalue
            maxkey = key
            maxvalue = value
        end
    end
    Pair(maxkey, maxvalue)
end
