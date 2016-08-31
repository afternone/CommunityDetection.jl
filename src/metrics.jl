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

function conductance(g, membership)
  nc = maximum(membership)
  ms = zeros(Int, nc)
  cs = zeros(Int, nc)
  ss = zeros(Int, nc)
  for u in vertices(g)
    u_idx = vertex_index(u,g)
    u_comm = membership[u_idx]
    ss[u_comm] += 1
    for v in out_neighbors(u,g)
      v_idx = vertex_index(v,g)
      v_comm = membership[v_idx]
      if u_comm == v_comm
        ms[u_comm] += 1
      else
        cs[u_comm] += 1
      end
    end
  end
  ss, [ms[i]+cs[i] > 0 ? cs[i]/(ms[i]+cs[i]) : 0. for i=1:nc]
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
