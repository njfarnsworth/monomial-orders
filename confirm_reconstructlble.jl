using LinearAlgebra

function matrix_validation(M)
    s = size(M);
    if s[1] != s[2] #checks that the matrix is square 
        return Inf
    end
    if Int(round(det(M))) == 0 # checks that the matrix is invertible
        return Inf
    end
    for col in 1:s[2] #checks that the first non-zero entry of every col > 0 
        for row in 1:s[1]
            if M[row, col] > 0
                break
            elseif M[row, col] < 0
                return Inf
            end
        end
    end
    return M
end 

function first_non_zero_entry(x) #returns the index and val of the first non-zero entry in a vector x
	for i in 1:length(x)
		if x[i] != 0
			return i, x[i]
		end
	end
	return Inf, Inf
end

function row_of_zeros(m) # iterates through an entire matrix and determines there there's a row of zeros somewhere
	for i in 1:size(m, 1)
		row = m[i,:]
		if first_non_zero_entry(row) == (Inf,Inf)
			return i
		end
	end
	return Inf
end

function ZARRE(M)
    n, _= size(M) # calculate how many rows the matrix has 

    for i in 1:n-1 # from the first row to the second to last row 
        M[i,:] = M[i,:] * (1/gcd(M[i,:]))
        col_index, _ = first_non_zero_entry(M[i,:]) # calculate the val/idx of the first non-zero entry in each row
        _, val = first_non_zero_entry(M[i,:])

        for j in i+1:n # move down one row, continue to the last row          
            if M[j, col_index] != 0 # check if the value under the first non-zero entry of the row is zero... 
                M[j,:] = abs(val) * M[j,:] - (val * M[j, col_index])/abs(val) * M[i,:]
            end
        end
    end

    return M 
end

function upper_triangularize(M) 
    n, _ = size(M)
    result = zeros(Int64,n,n)

    for i in 1:n
        first_non_zero_col, _ = first_non_zero_entry(M[i,:])
        result[:,i] = M[:,first_non_zero_col]
    end

    return result
end 

function reconstructible(M)
    n, _ = size(M)

    for i in 1:n-2
        A = M[:, 1:end .!= i]
        idx, _ = first_non_zero_entry(A[i,:])

        if idx == Inf # if the matrix now has a row of zeroes, return true 
            return true

        elseif n > 3
            c = A[idx+1, idx] / A[i, idx]

            for j in idx+1:n
                if A[idx+1, j] == c * A[i,j] #if multiple
                    continue
                else
                    return false
                end
            end

         else 
            return false
         end 
    end
    return true
end

M = [1 0 4; 0 1 0; 0 0 1] # Input Matrix Here 

if matrix_validation(M) != Inf
    print(reconstructible(upper_triangularize(ZARRE(M))))
end