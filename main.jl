using DataFrames, CSV

#TODO: make code to get a suitable equation of state given e and P or, alternatively, a CompOSE eos table
#TODO: eventually, remove the logs
#TODO: add all code mentions
#TODO: make include directory for the exported C++ code

begin

"""
    filter_decreasing(p, eps)

Filters two vectors, `p` and `eps`, to ensure that the `eps` vector
is monotonically increasing. It discards elements from both vectors
where a value in `eps` is smaller than the maximum `eps` value seen so far.

Returns a tuple containing the two filtered vectors: `(p_filtered, eps_filtered)`.
"""
function filter_decreasing(p::Vector, eps::Vector)
    # Return empty vectors if input is empty or lengths mismatch
    if isempty(eps) || length(p) != length(eps)
        return similar(p, 0), similar(eps, 0)
    end

    # Initialize the filtered vectors with the first element
    p_filtered = [p[1]]
    eps_filtered = [eps[1]]

    # Keep track of the last valid maximum value from eps
    last_max_eps = eps[1]

    # Iterate through the vectors starting from the second element
    for i in 2:length(eps)
        # Only keep the point if the eps value is not decreasing
        if eps[i] >= last_max_eps
            push!(p_filtered, p[i])
            push!(eps_filtered, eps[i])
            # Update the last maximum value
            last_max_eps = eps[i]
        end
        # If eps[i] < last_max_eps, the point is discarded
    end

    return p_filtered, eps_filtered
end
    
eosdf = CSV.read("testeos2.csv", DataFrame)
p = eosdf[:,1]
eps = eosdf[:,2]
p, eps = filter_decreasing(p, eps)
eps, p = filter_decreasing(eps, p)
@assert length(eps) == length(p)
@assert issorted(eps) && issorted(p)
n = length(eps)
out = zeros(2)

function getMR(eps, p, n, epsc)
    out = zeros(2)
    ccall((:qlimr_getMR, "./libqlimr.so"), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ptr{Cdouble}), p, eps, n, epsc, out)
    return copy(out)
end

end

epsc = exp.(range(log(10.0), log(eps[end]-100.0), length=100))
M = zeros(length(epsc))
R = zeros(length(epsc))
for i in 1:length(epsc)
    out = getMR(eps, p, n, epsc[i])
    M[i] = out[1]
    R[i] = out[2] 
end

df = DataFrame(epsc=epsc, M=M, R=R)
CSV.write("MRcpp2.csv", df, writeheader=false)