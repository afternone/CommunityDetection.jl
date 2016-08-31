"The binary Kullback-Leibler divergence."
function KL(q, p)
    (0.0 ≤ q ≤ 1.0 && 0.0 ≤ p ≤ 1.0) || error("q and p must in [0,1]")
    kl = 0.0
    if (q > 0.0 && p > 0.0)
        kl += q*log(q/p)
    end
    if q < 1.0 && p < 1.0
        kl += (1.0-q)*log((1.0-q)/(1.0-p))
    end
    kl
end

function permute_labels!(membership::Vector{Int})
    N = length(membership)
    if maximum(membership) > N || minimum(membership) < 1
        error("Label must between 1 and |V|")
    end
    label_counters = zeros(Int, N)
    j = 1
    for i=1:length(membership)
        k = membership[i]
        if k >= 1
            if label_counters[k] == 0
                # We have seen this label for the first time
                label_counters[k] = j
                k = j
                j += 1
            else
                k = label_counters[k]
            end
        end
        membership[i] = k
    end
end

"""get all groups"""
function getgrp(node_memory::Vector{Dict{Int,Int}})
    groups = Dict{Int,Vector{Int}}()
    for i=1:length(node_memory)
        for k in keys(node_memory[i])
            if haskey(groups, k)
                push!(groups[k], i)
            else
                groups[k] = [i]
            end
        end
    end
    collect(values(groups))
end

"""
calculate the entropy of a probility distribution
"""
function entropy(P::Vector{Float64})
    H = 0.0
    for i=1:length(P)
        # treat 0*log(0) as being equal to zero
        H += P[i] > 0.0 ? P[i]*log(P[i]) : 0.0
    end
    -H > 0 ? -H : 0.0
end

function entropy(membership::Vector{Int})
    n = length(membership)
    nb_comms = maximum(membership)
    groups = Array(Set{Int}, nb_comms)
    for i=1:nb_comms
        groups[i] = Set{Int}()
    end
    for i=1:n
        push!(groups[membership[i]], i)
    end
    h = 0.0
    for grp in groups
        n_c = length(grp)
        h += -plogp(n_c/n)
    end
    h
end

function entropy{V}(partition::AbstractPartition{V})
    g = graph(partition)
    n = num_vertices(g)
    h = 0.0
    for grp in values(partition.community)
        n_c = length(grp.nodes)
        h += -plogp(n_c/n)
    end
    h
end

"""read groups from file"""
function readgrp(filename)
    groups = Vector{Int}[]
    f = open(filename, "r")
    for ll in eachline(f)
        push!(groups, [parse(Int, i) for i in split(chomp(ll))])
    end
    groups
end

"""write groups to file"""
function writegrp(filename, groups::Vector{Vector{Int}})
    f = open(filename, "w")
    for i=1:length(groups)
        for j=1:length(groups[i])-1
            print(f, groups[i][j],' ')
        end
        print(f, groups[i][end],'\n')
    end
    close(f)
end

"""read membership from file"""
function readmsp(filename)
    membership = Dict{Int, Vector{Int}}()
    f = open(filename, "r")
    for ll in eachline(f)
        entries = [parse(Int, i) for i in split(chomp(ll))]
        membership[entries[1]] = entries[2:end]
    end
    membership
end

"""write membership to file"""
function writemsp(filename, membership::Dict{Int, Vector{Int}})
    f = open(filename, "w")
    for k in sort(collect(keys(membership)))
        print(f, k, '\t')
        for j=1:length(membership[k])-1
            print(f, membership[k][j], ' ')
        end
        print(f, membership[k][end], '\n')
    end
    close(f)
end

"""transform membership to groups"""
function msp2grp(membership::Dict{Int, Vector{Int}})
    groups = Dict{Int, Vector{Int}}()
    for (k,v) in membership
        for i in v
            if haskey(groups, i)
                push!(groups[i], k)
            else
                groups[i] = [k]
            end
        end
    end
    collect(values(groups))
end

"""transform groups to membership"""
function grp2msp(groups::Vector{Vector{Int}})
    membership = Dict{Int, Vector{Int}}()
    for i=1:length(groups)
        for j=1:length(groups[i])
            if haskey(membership, groups[i][j])
                push!(membership[groups[i][j]], i)
            else
                membership[groups[i][j]] = [i]
            end
        end
    end
    membership
end

# similarity between two nodes of an edge, (num_common_neighbors+1) / (num_total_neighbors+1)
function similarity1{V,E}(g::AbstractGraph{V,E}, e::E, neivec::Vector{Bool})
    u = source(e, g)
    v = target(e, g)
    for u_nei in out_neighbors(u, g)
        neivec[vertex_index(u_nei,g)] = true
    end
    num_common_nei = 0
    for v_nei in out_neighbors(v, g)
        if neivec[vertex_index(v_nei,g)]
            num_common_nei += 1
        end
    end
    # reset neivec
    for u_nei in out_neighbors(u, g)
        neivec[vertex_index(u_nei,g)] = false
    end

    num_total_nei = out_degree(u,g) + out_degree(v,g) - num_common_nei
    (num_common_nei + 2) / num_total_nei
end

function similarity1{V,E}(g::AbstractGraph{V,E})
    neivec = fill(false, num_vertices(g))
    Float64[similarity(g,e,neivec) for e in edges(g)]
end

function similarity{V,E}(g::AbstractGraph{V,E}, e::E, neivec::Vector{Bool})
    u = source(e, g)
    v = target(e, g)
    u_deg = out_degree(u,g)
    v_deg = out_degree(v,g)
    if u_deg > v_deg
    	u, v = v, u
    end
    for u_nei in out_neighbors(u, g)
        neivec[u_nei] = true
    end
    num_common_nei = 0
    for v_nei in out_neighbors(v, g)
        if neivec[v_nei]
            num_common_nei += 1
        end
    end
    # reset neivec
    for u_nei in out_neighbors(u, g)
        neivec[u_nei] = false
    end

    num_total_nei = u_deg + v_deg - num_common_nei
    (num_common_nei + 2) / num_total_nei
end

function similarity{V,E}(g::AbstractGraph{V,E})
    neivec = fill(false, num_vertices(g))
    Float64[similarity(g,e,neivec) for e in edges(g)]
end

# number of common neighbors between two nodes of an edge
function num_common_neighbors{V,E}(g::AbstractGraph{V,E}, e::E)
  u = source(e, g)
  v = target(e, g)
  u_nei = out_neighbors(u, g)
  v_nei = out_neighbors(v, g)
  length(intersect(u_nei, v_nei))
end

num_common_neighbors{V,E}(g::AbstractGraph{V,E}) = Int[num_common_neighbors(g,e) for e in edges(g)]
