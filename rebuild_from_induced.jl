include("confirm_reconstructlble.jl")


function induced_matrices_validation(matrices)
    for m in matrices
        rows, _ = size(m)
        if rows != length(matrices)-1 || matrix_validation(m) == Inf # n matrices of length n-1 x n-1 
            return Inf 
        end
    end

    return matrices
end

function gcd_matrix(m)
    m .= [(row / gcd(row)) for row in eachrow(m)]
    return m
end

function rebuildable(m1, m2)
    row_zeros_m1 = row_of_zeros(m1)
    row_zeros_m2 = row_of_zeros(m2)
    if row_zeros_m1 != Inf && row_zeros_m2 != Inf
        m1_no_row_zeros = m1[1:end .!= row_zeros_m1,:]
        m2_no_row_zeros = m2[1:end .!= row_zeros_m2, :]
        if gcd_matrix(m1_no_row_zeros) == gcd_matrix(m2_no_row_zeros)
            return true
        end
    end
    return false
end

function c_maker(matrix, n)
    c = zeros(Float64, 5, n) 

    for i in 1:n
        c[1, i] = i
        _, c[2, i] = first_non_zero_entry(matrix[:,i])
        
        c3_index, _ = first_non_zero_entry(matrix[c[2,i],:])
        c[3, i] = matrix[c[2,i], c3_index] 

        idx_ignore, _ = first_non_zero_entry(matrix[c[2,i], :])
        _, c[4, i] = first_non_zero_entry(matrix[c[2,i], idx_ignore+1:end])
        c[5, i], _ = first_non_zero_entry(matrix[c[2,i], idx_ignore+1:end])
    end
    return c
end

function row_operation(row1, row2)
    new_row = []

    _, row1_non_zero_val = first_non_zero_entry(row1)
    _, row2_non_zero_val = first_non_zero_entry(row2)

    if row1_non_zero_val * row2_non_zero_val > 0
        new_row = row1 - row2
    else 
        new_row = row1 + row2
    end
    return new_row
end

function induced_matrices(matrix, n) 
    induced_matrices_list = Matrix{Float64}[]
    C = c_maker(matrix, n) 
    for i in 1:n
        A = matrix[:, 1:end .!= i ]

        s = C[2, i]
        a = C[4, i]

        if a!= 0
            e = C[2, C[5, i]]
            b = C[3,C[5,i]]
        end

        while  a!= 0  
            A[e, :] = row_operation(A[e,:]* abs(a), (A[s,:] * abs(b)) * (a*b)/abs(a*b))
            s = e
            _, a = first_non_zero_entry(A[s,:]) 
            c, _ = first_non_zero_entry(A[s,:])
            e = C[2, c+1]
            b = C[3, c+1]     
        end

        final = A[1:end .!= s, :]
        
        for i in 1:n-1
            d = gcd(final[i,:])
            final[i,:] = final[i,:] * (1/d)
        end

        push!(induced_matrices_list, final)
    end

    return induced_matrices_list

end

function reconstruction(matrices)
    n = size(matrices[1])[1] + 1 # size of the reconstructed matrix
    m = zeros(Int, n,n)
  
    second_to_last_m_last_col_removed = matrices[n-1][:, 1:end .!= n-1]
    last_m_last_col_removed = matrices[n][:, 1:end .!= n-1]

    if rebuildable(second_to_last_m_last_col_removed, last_m_last_col_removed)
        print("There is at least one MORMORE corresponding to these induced orders")
         # 4 cases total
        rozsl = row_of_zeros(second_to_last_m_last_col_removed)
        rozl = row_of_zeros(last_m_last_col_removed)

        if rozsl < rozl # case 1
            m[rozsl, n] = matrices[n-1][rozsl, n-1] 
            m[rozl+1, n-1] = matrices[n][rozl, n-1]

            for i = in 1:n

                if i < rozsl
                    _, first_non_zero_val_stl = first_non_zero_entry(matrices[n-1][i,:])
                    _, first_non_zero_val_last = first_non_zero_entry(matrices[n][i,:])

                    for j in 1:n-1
                        m[i,j] = (abs(first_non_zero_val_stl)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n][i,j]
                    end
                    m[i,n] = (abs(first_non_zero_val_last)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n-1][i,n-1]

                elseif i > rozsl && i <= rozl
                    _, first_non_zero_val_stl = first_non_zero_entry(matrices[n-1][i,:])
                    _, first_non_zero_val_last = first_non_zero_entry(matrices[n][i-1,:])

                    for j in 1:n-1
                        m[i,j] = (abs(first_non_zero_val_stl)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n][i-1,j]
                    end
                    m[i,n] = (abs(first_non_zero_val_last)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n-1][i,n-1]

                elseif i > rozl+1
                    _, first_non_zero_val_stl = first_non_zero_entry(matrices[n-1][i-1,:])
                    _, first_non_zero_val_last = first_non_zero_entry(matrices[n][i-1,:])

                    for j in 1:n-1
                        m[i,j] = (abs(first_non_zero_val_stl)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n][i-1,j]
                    end
                    m[i,n] = (abs(first_non_zero_val_last)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n-1][i-1,n-1]  
                end
            end
        end

        # case 1 complete 

        elseif rozsl > rozl # case 2 
            m[rozsl+1, n] = matrices[n-1][rozsl, n-1] 
            m[rozl,n-1] = matrices[n][rozl, n-1]

            for i in 1:n
                if i < rozl
                    _, first_non_zero_val_stl = first_non_zero_entry(matrices[n-1][i,:])
                    _, first_non_zero_val_last = first_non_zero_entry(matrices[n][i,:])

                    for j in 1:n-1
                        m[i,j] = (abs(first_non_zero_val_stl)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n][i,j]
                    end
                    m[i,n] = (abs(first_non_zero_val_last)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n-1][i,n-1]
                
                elseif i > rozl && i < rozsl + 1
                    _, first_non_zero_val_stl = first_non_zero_entry(matrices[n-1][i-1,:])
                    _, first_non_zero_val_last = first_non_zero_entry(matrices[n][i,:])

                    for j in 1:n-1
                        m[i,j] = (abs(first_non_zero_val_stl)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n][i,j]
                    end
                    m[i,n] = abs(first_non_zero_val)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last)) * matrices[n-1][i-1,n-1]

                elseif i > rozsl + 1
                    _, first_non_zero_val_stl = first_non_zero_entry(matrices[n-1][i-1,:])
                    _, first_non_zero_val_last = first_non_zero_entry(matrices[n][i-1,:])
            
                    for j in 1:n-1
                        m[i,j] = (abs(first_non_zero_val_stl)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last))) * matrices[n][i-1,j]
                    end
                    m[i,n] = abs(first_non_zero_val)/gcd(abs(first_non_zero_val_stl),abs(first_non_zero_val_last)) * matrices[n-1][i-1,n-1]
                end
            end
   
        elseif rozsl == rozl && rozl != n-1 # case 3
            s, _ = first_non_zero_entry(matrices[n][rozl+1, :])
            s_induced = matrices[s]
            m[rozl, n-1] = s_induced[rozl,n-2]
            m[rozl, n] = s_induced[rozl, n-1]

            for i in 1:rozl-1
                _, fnzl = first_non_zero_entry(matrices[n][i,:])
                _, fnzsl = first_non_zero_entry(matrices[n-1][i,:])

                for j in 1:n-1
                    m[i,j] = abs(fnzsl)/gcd(abs(fnzl),abs(fnzsl)) * matrices[n][i,j]
                end
                m[i,n] = abs(fnzl)/gcd(abs(fnzl),abs(fnzsl)) * matrices[n-1][i,n-1]
            end

            k = rozl + 1
            match = false  

            while (k <= n && !match)

                if s_induced[rozl,n-2] != 0
                    m[k,n] = 1             
                else
                    m[k,n-1] = 1       
                end 

                for i in rozl+1:k-1
                    _, fnzl = first_non_zero_entry(matrices[n][i-1,:])
                    _, fnzsl = nonZeroVal(matrices[n-1][i,:])

                    for j in 1:n-1
                        m[i,j] = abs(fnzsl)/gcd(abs(fnzl),abs(fnzsl)) * matrices[n][i-1,j]
                    end
                    m[i,n] = abs(fnzl)/gcd(abs(fnzl),abs(fnzsl)) * matrices[n-1][i,n-1]
                end

                for i in k+1:n
                    _, fnzl = first_non_zero_entry(matrices[n][i-1,:])
                    _, fnzsl = first_non_zero_entry(matrices[n-1][i-1,:])

                    for j in 1:n-1
                        m[i,j] = abs(fnzsl)/gcd(abs(fnzl),abs(fnzsl)) * matrices[n][i-1,j]
                    end
                    m[i,n] = abs(fnzl)/gcd(abs(fnzl),abs(fnzsl)) * matrices[n-1][i-1,n-1]


                end

                induced_matrices = induced_matrices(m,n) 
                
                if induced_matrices == matrices
                    match = true
                
                else
                    m[k:end,:] = zeros(Int,n-k,n)
                    k += 1
                end
            end  

        elseif rozsl == rozl && rozl == n-1 # case 4
            print("Unable to reconstruct matrix")
            return Inf

        
    else
        print("There is no MORMORE corresponding to these induced orders")
        return Inf
    end

    return m
    
end

# ---------------------------------------------------------------------------------------------

m1 = [1 4 5; 6 5 8; 3 8 9] #INPUT MATRICES HERE
m2 = [8 2 5; 3 4 5; 3 8 3]
m3 = [2 3 4; 5 7 5; 10 2 0]
m4 = [2 1 1; 5 6 7; 1 3 4]

matrices = [m1, m2, m3, m4] # MANUALLY GENERATE LIST 
matrices .= ZARRE.(matrices) 
result = reconstruction(matrices)


