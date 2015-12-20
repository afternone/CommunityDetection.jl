using CommunityDetection
using Base.Test
using Graphs

g = simple_graph(6, is_directed=false)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 1)
add_edge!(g, 1, 4)
add_edge!(g, 4, 5)
add_edge!(g, 5, 6)
add_edge!(g, 6, 4)

mp = mpartition(g)
optimize_partition(mp)
@test nmi(mp.membership, [1,1,1,2,2,2]) == 1.0
