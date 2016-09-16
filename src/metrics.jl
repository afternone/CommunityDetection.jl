plog2p(x) = x > 0.0 ? x*log2(x) : 0.0
function codelen(g,m)
  ne = num_edges(g)
  nc = maximum(m)
  exit_ne = zeros(Int,nc)
  inner_ne = zeros(Int,nc)
  total_exit_ne = 0
  for e in edges(g)
    u = source(e, g)
    v = target(e, g)
    u_idx = vertex_index(u, g)
    v_idx = vertex_index(v, g)
    u_comm = m[u_idx]
    v_comm = m[v_idx]
    if u_comm != v_comm
      exit_ne[u_comm] += 1
      exit_ne[v_comm] += 1
      total_exit_ne += 2
    else
      inner_ne[u_comm] += 2
    end
  end
  L1 = plog2p(total_exit_ne/2ne)
  L2 = -2sum([plog2p(x/2ne) for x in exit_ne])
  L3 = -sum([plog2p(out_degree(v,g)/2ne) for v in vertices(g)])
  L4 = sum([plog2p((exit_ne[i]+inner_ne[i])/2ne) for i=1:nc])
  L1 + L2 + L3 + L4
end

function min_conductance(g, m)
  x, y = conductance(g, m)
  nc_co = Dict{Int,Float64}()
  for i=1:length(x)
    if y[i] > 0.
      if haskey(nc_co, x[i])
        if nc_co[x[i]] > y[i]
          nc_co[x[i]] = y[i]
        end
      else
        nc_co[x[i]] = y[i]
      end
    end
  end
  xx = collect(keys(nc_co))
  yy = collect(values(nc_co))
  idx = sortperm(xx)
  xx[idx], yy[idx]
end

function min_cond(x,y)
  nc_co = Dict{Int,Float64}()
  for i=1:length(x)
    if y[i] > 0.
      if haskey(nc_co, x[i])
        if nc_co[x[i]] > y[i]
          nc_co[x[i]] = y[i]
        end
      else
        nc_co[x[i]] = y[i]
      end
    end
  end
  xx = collect(keys(nc_co))
  yy = collect(values(nc_co))
  idx = sortperm(xx)
  xx[idx], yy[idx]
end

"""
Get the number of internal nodes `ns`,
the number of internal edges `ms`,
and the number of external (boundary) edges `cs`.
"""
function ns_ms_cs(g, membership)
  nc = maximum(membership)
  ns = zeros(Int, nc)
  ms = zeros(Int, nc)
  cs = zeros(Int, nc)
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    ns[u_comm] += 1
    for v in out_neighbors(u,g)
      v_idx = vertex_index(v,g)
      v_comm = membership[v_idx]
      # count each edge once
      if u_idx < v_idx
        if u_comm == v_comm
          ms[u_comm] += 1
        else
          cs[u_comm] += 1
          cs[v_comm] += 1
        end
      end
    end
  end
  ns, ms, cs
end

function internal_density(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, [ns[i] > 1 ? 2*ms[i]/ns[i]/(ns[i]-1) : 0. for i=1:length(ns)]
end

function edges_inside(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, ms
end

function average_degree(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, 2ms./ns
end

""" Fraction over median degree is the number of nodes that
have an internal degree greater than the median degree of all nodes in the graph.
"""
function fomd(g, membership)
  median_deg = median([out_degree(u,g) for u in vertices(g)])
  nc = maximum(membership)
  ms = zeros(Int, nc)
  ns = zeros(Int, nc)
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    ns[u_comm] += 1
    i = 0
    for v in out_neighbors(u,g)
      v_idx = vertex_index(v,g)
      v_comm = membership[v_idx]
      if u_comm == v_comm
        i += 1
      end
    end
    if i > median_deg
      ms[u_comm] += 1
    end
  end
  ns, ms./ns
end

function clustering_coefficient(g, membership)
  nc = maximum(membership)
  ns = zeros(Int, nc)
  ms = zeros(Int, nc)
  cs = zeros(Int, nc)
  rv = fill(false, num_vertices(g))
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    ns[u_comm] += 1
    for v in out_neighbors(u,g)
      v_idx = vertex_index(v,g)
      v_comm = membership[v_idx]
      if u_idx < v_idx
        for w in out_neighbors(v,g)
          w_idx = vertex_index(w,g)
          w_comm = membership[w_idx]
          if v_idx < w_idx && u_comm == v_comm == w_comm
            if u in out_neighbors(w,g)
              ms[u_comm] += 1
            else
              cs[u_comm] += 2
            end
          end
        end
      end
    end
  end
  ns, [ms[i]/(ms[i]+cs[i]) for i=1:nc]
end

function triangle_participation(g)
  rv = fill(false, num_vertices(g))
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    rv[u_idx] && continue
    for v in out_neighbors(u,g)
      v_idx = vertex_index(v,g)
      for w in out_neighbors(v,g)
        w_idx = vertex_index(w,g)
        if u in out_neighbors(w,g)
          rv[u_idx] = true
          rv[v_idx] = true
          rv[w_idx] = true
        end
      end
    end
  end
  rv
end

function triangle_participation_ratio(g)
  rv = triangle_participation(g)
  sum(rv)/num_vertices(g)
end

function triangle_participation_ratio(g, membership)
  rv = triangle_participation(g)
  nc = maximum(membership)
  ms = zeros(Int, nc)
  ns = zeros(Int, nc)
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    ns[u_comm] += 1
    if rv[u_idx]
      ms[u_comm] += 1
    end
  end
  ns, ms./ns
end

function expansion(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, cs./ns
end

"""
Cut ratio is the ratio between the number of external (boundary) edges
in a cluster and the cluster's maximum possible number of external edges
"""
function cut_ratio(g, membership)
  n = num_vertices(g)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, [ns[i] < n ? cs[i]/(ns[i]*(n-ns[i])) : 0. for i=1:length(ns)]
end

"""
Conductance is the ratio between the number of external (boundary) edges
in a cluster and the cluster's total number of edges
"""
function conductance(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, [2ms[i]+cs[i] > 0 ? cs[i]/(2ms[i]+cs[i]) : 0. for i=1:length(ns)]
end

function normalized_cut(g, membership)
  m = num_edges(g)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, [2ms[i]+cs[i] > 0 && 2*(m-ms[i])+cs[i] > 0 ? cs[i]/(2ms[i]+cs[i])+cs[i]/(2*(m-ms[i])+cs[i]) : 0. for i=1:length(ns)]
end

""" out degree fraction """
function odf(g, membership)
  nc = maximum(membership)
  ns = zeros(Int, nc)
  degs = [out_degree(u,g) for u in vertices(g)]
  out_degs = zeros(Int, num_vertices(g))
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    ns[u_comm] += 1
    for v in out_neighbors(u,g)
      v_idx = vertex_index(v,g)
      v_comm = membership[v_idx]
      if u_idx < v_idx && u_comm != v_comm
        out_degs[u_idx] += 1
        out_degs[v_idx] += 1
      end
    end
  end
  ns, out_degs
end

""" maximum out degree fraction """
function max_odf(g, membership)
  ns, out_degs = odf(g, membership)
  out_deg_ratio = [out_degs[vertex_index(u,g)]/out_degree(u,g) for u in vertices(g)]
  cs = zeros(length(ns))
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    if cs[u_comm] < out_deg_ratio[u_idx]
      cs[u_comm] = out_deg_ratio[u_idx]
    end
  end
  ns, cs
end

function average_odf(g, membership)
  ns, out_degs = odf(g, membership)
  out_deg_ratio = [out_degs[vertex_index(u,g)]/out_degree(u,g) for u in vertices(g)]
  cs = zeros(length(ns))
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    cs[u_comm] += out_deg_ratio[u_idx]
  end
  ns, cs./ns
end

"""
flake_odf is the fraction of nodes in a community that have edges pointing
inside than to the outside of the cluster.
"""
function flake_odf(g, membership)
  ns, out_degs = odf(g, membership)
  cs = zeros(Int, length(ns))
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    if out_degree(u,g) < 2out_degs[u_idx]
      cs[u_comm] += 1
    end
  end
  ns, cs./ns
end

""" seprability is the ratio between the internal and the external number of edges """
function separability(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, ms./cs
end

"""
Density is the fraction of the edges (out of all possible edges)
that appear between the nodes in a community.
"""
function density(g, membership)
  ns, ms, cs = ns_ms_cs(g, membership)
  ns, [ns[i] > 1 ? 2ms[i]/ns[i]/(ns[i]-1) : 0. for i=1:length(ns)]
end

function cohesiveness(g, membership)
  # TO DO
end

function p_in_after_n(g, v, n, comm)
  p_in_after_n_r_cached(g, v, n, Set(comm), Dict{Pair{Int,Int},Float64}())
end

function p_in_after_n_r_cached(g, v, n, comm, cache)
  if haskey(cache, Pair(v,n))
    return cache[Pair(v,n)]
  end
  if !(v in comm)
    return 0.
  end
  neighbors = Set(out_neighbors(v,g))
  numNei = length(neighbors)
  if n == 1
    len_intersec = length(intersect(neighbors, comm))
    return len_intersec > 0 ? len_intersec/numNei : 0.
  end
  totalP = 0.
  for neighbor in neighbors
    pGivenNei = p_in_after_n_r_cached(g, neighbor, n-1, comm, cache)
    cache[Pair(neighbor, n-1)] = pGivenNei
    totalP += 1/numNei*pGivenNei
    cache[Pair(v,n)] = totalP
    return totalP
  end
  0.
end

function prob_metric_cluster(g, members)
  nc = length(members)
  data = [p_in_after_n(g,v,nc,members) for v in members]
  mean_ = mean(data)
  std_ = std(data)
  var_ = var(data)
  mean_, std_, var_
end

"""
    Calculates the probability metric on the graph G for each cluster in
    clusters. Returns a list of 3-tuples [(a, b, c),...] where a is the mean,
    b the standard deviation, and c the variance, indexed by cluster id.
    This metric measures how likely a particle placed on some vertex will stay within
    the original community after n random steps, where n is the number of vertices in
    the community (or some other, better value for normalization).
    This returns a list [(mean, variance)_0, (mean, variance)_1..., (mean, variance)_n]
    where the statistics within each cluster are represented by its index in the list.
    A high mean suggests strong community structure, but a high variance suggests that
    a few objects in the community might be outliers and not necessarily belong.

    Currently, we use the size of the community as the number of steps taken, but this
    isn't based on much theory. Some analysis is necessary to find a better choice for n.
    If you want one "score" to score the clustering instead of individual communities,
    call probability_metric_score instead.
"""
function prob_metric_graph(g, membership)
  nc = maximum(membership)
  full_mean = zeros(nc)
  full_std = zeros(nc)
  full_var = zeros(nc)
  members = [Vector{Int}() for i=1:nc]
  for i=1:length(membership)
    push!(members[membership[i]], i)
  end
  for i=1:nc
    full_mean[i], full_std[i], full_var[i] = prob_metric_cluster(g, members[i])
  end
  full_mean, full_std, full_var
end
