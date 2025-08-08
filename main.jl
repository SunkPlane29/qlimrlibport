using DataFrames, CSV

#TODO: make code to get a suitable equation of state given e and P or, alternatively, a CompOSE eos table
#TODO: eventually, remove the logs
#TODO: add all code mentions
#TODO: make include directory for the exported C++ code

begin

"""
    filter_both_sorted(p, eps)

Filters two vectors, `p` and `eps`, keeping only the data points (p[i], eps[i])
that maintain a monotonically increasing order in BOTH vectors simultaneously.

Returns a tuple containing the two filtered vectors: `(p_filtered, eps_filtered)`.
"""
function filter_both_sorted(p::Vector, eps::Vector)
    # Return empty vectors if input is empty or lengths mismatch
    if isempty(p) || length(p) != length(eps)
        return similar(p, 0), similar(eps, 0)
    end

    # Initialize the filtered vectors with the first element
    p_filtered = [p[1]]
    eps_filtered = [eps[1]]

    # Keep track of the last maximum value seen in BOTH vectors
    last_max_p = p[1]
    last_max_eps = eps[1]

    # Iterate through the vectors starting from the second element
    for i in 2:length(p)
        # The condition now checks both vectors
        if p[i] >= last_max_p && eps[i] >= last_max_eps
            push!(p_filtered, p[i])
            push!(eps_filtered, eps[i])

            # Update both "high-water marks"
            last_max_p = p[i]
            last_max_eps = eps[i]
        end
        # If the condition is false, the point is discarded
    end

    return p_filtered, eps_filtered
end
    
eosdf = CSV.read("testeos2.csv", DataFrame)
p = eosdf[:,1]
eps = eosdf[:,2]
p, eps = filter_both_sorted(p, eps)
@assert length(eps) == length(p)
@assert issorted(eps) && issorted(p)
n = length(eps)

function getMR(eps, p, n, epsc, R_initial)
    out = zeros(2)
    ccall((:qlimr_getMR, "./libqlimr.so"), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Cdouble, Ptr{Cdouble}), p, eps, n, epsc, R_initial, out)
    return copy(out)
end

end

# getMR(eps, p, n, 300.0, 0.0004)

epsc = exp.(range(log(10.0), log(eps[end]-100.0), length=100))
M = zeros(length(epsc))
R = zeros(length(epsc))
@time for i in 1:length(epsc)
    out = getMR(eps, p, n, epsc[i], 0.0004)
    M[i] = out[1]
    R[i] = out[2] 
end

df = DataFrame(epsc=epsc, M=M, R=R)
CSV.write("MRcpp2.csv", df, writeheader=false)