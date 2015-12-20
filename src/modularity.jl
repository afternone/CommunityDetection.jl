modularity{V}(graph::AbstractGraph{V}, membership::Vector{Int}) = modularity(graph, membership, Float64[])

function modularity{V,T<:Real}(graph::AbstractGraph{V}, membership::Vector{Int}, weights::Vector{T})

    types = maximum(membership)
    no_of_edges = num_edges(graph)
    m = 0.0

    if length(membership) < num_vertices(graph)
        error("cannot calculate modularity, membership vector too short")
    end

    e = zeros(types)
    a = zeros(types)

    if !isempty(weights)
        if length(weights) < no_of_edges
            error("cannot caculate modularity, weight vector too short")
        end
        m = sum(weights)
        for i in edges(graph)
            w = weights[edge_index(i, graph)]
            w >= 0 || error("negative weight in weight vector")
            c1 = membership[vertex_index(source(graph, i), graph)]
            c2 = membership[vertex_index(target(graph, i), graph)]
            if c1 == c2
                e[c1] += 2w
            end
            a[c1] += w
            a[c2] += w
        end
    else
        m = no_of_edges
        for i in edges(graph)
            c1 = membership[vertex_index(source(i, graph), graph)]
            c2 = membership[vertex_index(target(i, graph), graph)]
            if c1 == c2
                e[c1] += 2
            end
            a[c1] += 1
            a[c2] += 1
        end
    end

    modularity_value = 0.0
    if m > 0
        for i=1:types
            tmp = a[i]/2/m
            modularity_value += e[i]/2/m
            modularity_value -= tmp*tmp
        end
    end
    modularity_value
end
