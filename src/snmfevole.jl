function snmf_evolve!(A, K, δ, T, D, H, α)
    n = size(A,1)
    B = H*D
    for t=1:T
        Z = H*D*H'
        logl = α*dKL_mat(A, Z) + (1-α)*dKL_mat(B, H*D)
        for k=1:K, i=1:n
            s = 0.
            for j=1:n
                if Z[i,j] > 0
                    s += A[i,j]*D[k,k]*H[j,k]/Z[i,j]
                end
            end
            H[i,k] = H[i,k]*2α*s + (1-α)*B[i,k]
        end
        H ./= sum(H,1)

        Z = H*D*H'
        for k=1:K
            s1 = 0.
            s2 = 0.
            for i=1:n
                s1 += B[i,k]
                for j=1:n
                    if Z[i,j] > 0
                        s2 += A[i,j]*H[i,k]*H[j,k]/Z[i,j]
                    end
                end
            end
            D[k,k] = D[k,k]*α*s2 + (1-α)*s1
        end
        D ./= sum(diag(D))

        logl1 = α*dKL_mat(A, H*D*H') + (1-α)*dKL_mat(B, H*D)
        if  abs(logl1 - logl) <= δ
            return logl1
        end
    end
end

function snmf_evolve(A, K, δ, T)
    n = size(A,1)
    H = rand(n,K)
    D = diagm(rand(K))
    H ./= sum(H,1)
    D ./= sum(diag(D))
    A ./= sum(A)
    loglike = dKL_mat(A, H*D*H')
    for t=1:T
        Z = H*D*H'
        logl = dKL_mat(A, Z)
        for k=1:K, i=1:n
            s = 0.
            for j=1:n
                if Z[i,j] > 0
                    s += A[i,j]*D[k,k]*H[j,k]/Z[i,j]
                end
            end
            H[i,k] = H[i,k]*2*s
        end
        H ./= sum(H,1)

        Z = H*D*H'
        for k=1:K
            s = 0.
            for i=1:n, j=1:n
                if Z[i,j] > 0
                    s += A[i,j]*H[i,k]*H[j,k]/Z[i,j]
                end
            end
            D[k,k] = D[k,k]*s
        end
        D ./= sum(diag(D))

        logl1 = dKL_mat(A, H*D*H')
        println(t, '\t', abs(logl1-logl))
        if  abs(logl1 - logl) <= δ
            loglike = logl1
            break
        end
    end
    D, H, loglike
end

function dKL_mat(X,Y)
    X=X./sum(X)
    Y=Y./sum(Y)
    X = X+eps()
    Y = Y+eps()
    sum(X .* log(X ./ Y))
end
