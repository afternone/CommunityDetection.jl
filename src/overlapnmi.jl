"""
Normalized mutual information measure
Measures the similarity of two covers, output a float number in range [0, 1],
with 1 corresponding to a perfect matching.
c1 => partition 1, c2 => partition 2, N => number of total nodes
"""
function nmi(c1::Vector{Vector{Int}}, c2::Vector{Vector{Int}}, N::Int=maximum(map(maximum, c1)))
    k, l = length(c1), length(c2)
    Hx_vector = zeros(k) # vector with length k, entrophy of each cluster in clusters1
    Hy_vector = zeros(l) # vector with length l, entrophy of each cluster in clusters2
    Hxy_martix = zeros(k, l) # 2D array, k*l, relative entropy H(x_k,y_l)
    Hyx_martix = zeros(l, k) # 2D array, l*k, relative entropy H(y_l,x_k)
    for i=1:k
        px1 = length(c1[i])/N # probability of a random node belongs to this cluster
        Hx_vector[i] = entropy([px1,1.0-px1])
    end
    for j=1:l
        py1 = length(c2[j])/N
        Hy_vector[j] = entropy([py1,1.0-py1])
    end
    for i=1:k, j=1:l
        px1y1 = length(intersect(c1[i], c2[j]))/N # joint probability
        px1y0 = length(setdiff(c1[i], c2[j]))/N
        px0y1 = length(setdiff(c2[j], c1[i]))/N
        px0y0 = 1.0 - px1y1 - px1y0 - px0y1
        joint_Hxy = entropy([px1y1, px1y0, px0y1, px0y0])
        Hxy_martix[i,j] = joint_Hxy - Hy_vector[j]
        Hyx_martix[j,i] = joint_Hxy - Hx_vector[i]
    end
    Hxy_vector = collect(minimum(Hxy_martix, 2))
    Hyx_vector = collect(minimum(Hyx_martix, 2))
    Hxy = 0.0
    Hyx = 0.0
    for i=1:k
        Hxy += Hx_vector[i] > 0.0 ? Hxy_vector[i]/Hx_vector[i] : 0.0
    end
    Hxy /= k
    for j=1:l
        Hyx += Hy_vector[j] > 0.0 ? Hyx_vector[j]/Hy_vector[j] : 0.0
    end
    Hyx /= l
    1.0 - 0.5*(Hxy + Hyx)
end
