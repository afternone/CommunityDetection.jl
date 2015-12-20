"""
This function finds the most frequent labels of neighbors
and randomly returns one of them
"""
function maxvote{T}(label_count::Dict{Int,T})
    max_count = maximum(values(label_count))
    candidate_labels = Int[]
    for (k,v) in label_count
        if v == max_count
            push!(candidate_labels, k)
        end
    end
    sample(candidate_labels)
    #candidate_labels[rand(1:length(candidate_labels))]
end

"""
Applies the value of label for each neighbor
and calls maxVote function
Use multinomial sampling for speaker rule
Use maximum vote for listener rule
"""
function applyvote!{V}(v::V, g::AbstractGraph{V}, node_memory::Vector{Dict{Int,Int}}, β::Float64=1.0)
    v_neighbors = out_neighbors(v, g)
    v_idx = vertex_index(v, g)
    if length(v_neighbors)>0
        label_list = Dict{Int, Int}()
        for u in v_neighbors
            u_idx = vertex_index(u, g)
            label = sample(collect(keys(node_memory[u_idx])), WeightVec(collect(values(node_memory[u_idx])).^β))
            if haskey(label_list, label)
                label_list[label] += 1
            else
                label_list[label] = 1
            end
        end
        # listener chose a received label to add to memory
        selected_label = maxvote(label_list)
        # add the selected label to the memory
        if haskey(node_memory[v_idx], selected_label)
            node_memory[v_idx][selected_label] += 1
        else
            node_memory[v_idx][selected_label] = 1
        end
    end
end

"""
Applies the value of label for each neighbor
and calls maxVote function
Use multinomial sampling for speaker rule
Use neighbor strength driven for listener rule
"""
function applyvote!{V}(v::V, g::AbstractGraph{V}, node_memory::Vector{Dict{Int,Int}}, ns::Vector{Int}, β::Float64=1.0)
    v_edges = out_edges(v, g)
    v_idx = vertex_index(v, g)
    if length(v_edges)>0
        label_list = Dict{Int, Int}()
        for e in v_edges
            e_idx = edge_index(e, g)
            u = target(e, g)
            u_idx = vertex_index(u, g)
            label = sample(collect(keys(node_memory[u_idx])), WeightVec(collect(values(node_memory[u_idx])).^β))
            if haskey(label_list, label)
                label_list[label] += 1 + ns[e_idx]
            else
                label_list[label] = 1 + ns[e_idx]
            end
        end
        # listener chose a received label to add to memory
        selected_label = maxvote(label_list)
        # add the selected label to the memory
        if haskey(node_memory[v_idx], selected_label)
            node_memory[v_idx][selected_label] += 1
        else
            node_memory[v_idx][selected_label] = 1
        end
    end
end

function slpa{V}(g::AbstractGraph{V}, β::Float64=1.0, iterations::Int=20)
    node_memory = Dict{Int,Int}[]
    for i=1:num_vertices(g)
        push!(node_memory, Dict{Int,Int}(i=>1))
    end
    for i=1:iterations
        order = shuffle(collect(vertices(g)))
        for v in order
            applyvote!(v, g, node_memory, β)
        end
    end
    node_memory
end

function slpa{V}(g::AbstractGraph{V}, ns::Vector{Int}, β::Float64=1.0, iterations::Int=20)
    node_memory = Dict{Int,Int}[]
    for i=1:num_vertices(g)
        push!(node_memory, Dict{Int,Int}(i=>1))
    end
    for i=1:iterations
        # shuffe nodes order
        order = shuffle(collect(vertices(g)))
        for v in order
            applyvote!(v, g, node_memory, ns, β)
        end
    end
    node_memory
end

function slpa{V}(g::AbstractGraph{V}; r::Float64=0.2, β::Float64=1.0, iterations::Int=20)
    node_memory = slpa(g, β, iterations)
    post_processing!(node_memory, r)
    node_memory
end

function slpa{V}(g::AbstractGraph{V}, ns::Vector{Int}; r::Float64=0.2, β::Float64=1.0, iterations::Int=20)
    node_memory = slpa(g, ns, β, iterations)
    post_processing!(node_memory, r)
    node_memory
end

# performs post processing to remove the labels that are below the threshhold
# r∈[0,1], if the probability is less than r, remove it during post processing
function post_processing!(node_memory::Vector{Dict{Int,Int}}, r::Real)
    for memory in node_memory
        tempmemory = copy(memory)
        sum_count = sum(values(memory))
        threshold = sum_count*r
        for (k,v) in memory
            if v < threshold
                delete!(memory, k)
            end
        end

        # if r is too large to lead to some node_memory empty, we will select label with the max_count
        if isempty(memory)
            selected_label = maxvote(tempmemory)
            memory[selected_label] = tempmemory[selected_label]
        end
    end
end
