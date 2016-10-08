function modularity_density(g,c)
	m = num_edges(g)
	nc = maximum(c)
	ns = zeros(Int,nc)
	Ein = zeros(Int,nc)
	Eout = zeros(Int,nc)
	g_comm = [Dict{Int,Int}() for i=1:nc]
	for u in vertices(g)
		u_idx = vertex_index(u,g)
		u_comm = c[u_idx]
		ns[u_comm] += 1
	end
	for edge in edges(g)
		u,v = source(edge,g),target(edge,g)
		u_idx,v_idx = vertex_index(u,g),vertex_index(v,g)
		u_comm,v_comm = c[u_idx],c[v_idx]
		if u_comm == v_comm
			Ein[u_comm] += 1
		else
			Eout[u_comm] += 1
			Eout[v_comm] += 1
			if haskey(g_comm[u_comm],v_comm)
				g_comm[u_comm][v_comm] += 1
			else
				g_comm[u_comm][v_comm] = 1
			end
			if haskey(g_comm[v_comm],u_comm)
				g_comm[v_comm][u_comm] += 1
			else
				g_comm[v_comm][u_comm] = 1
			end
		end
	end
	Qds = 0.
	for i=1:nc
		dci = ns[i] > 1 ? 2Ein[i]/ns[i]/(ns[i]-1) : 0.
		Qds += dci*Ein[i]/m - (dci*(2Ein[i]+Eout[i])/2m)^2
		for nb_comm in keys(g_comm[i])
			Qds -= g_comm[i][nb_comm]^2/2m/ns[i]/ns[nb_comm]
		end
	end
	Qds
end