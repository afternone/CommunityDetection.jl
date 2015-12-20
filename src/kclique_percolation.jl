function kclique_percolation(g::SimpleGraph, k::Int)
    kclqs = filter(x->length(x)>=k, maximal_cliques(g))
    n_c = length(kclqs)
    h = Graph(n_c)
    
    for i=1:n_c, j=i+1:n_c
        if length(intersect(kclqs[i], kclqs[j]))>=k-1
            add_edge!(h, i, j)
        end
    end
    
    cps = connected_components(h)
    comms = Array(Set{Int}, length(cps))
    for i=1:length(comms)
        comms[i] = Set{Int}()
    end
    
    for i=1:length(cps)
        for cp in cps[i]
            push!(comms[i], kclqs[cp]...)
        end
    end
    comms
end