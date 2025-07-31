using DataFrames, CSV

#TODO: make code to get a suitable equation of state given e and P or, alternatively, a CompOSE eos table
#TODO: eventually, remove the logs
#TODO: add all code mentions
#TODO: make include directory for the exported C++ code

begin
    
eosdf = CSV.read("testeos.csv", DataFrame)
eps = eosdf[:,1]
p = eosdf[:,2]
@assert length(eps) == length(p)
n = length(eps)
out = zeros(2)

function getMR(eps, p, n, epsc)
    out = zeros(2)
    ccall((:qlimr_getMR, "./libqlimr.so"), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ptr{Cdouble}), p, eps, n, epsc, out)
    return copy(out)
end

end

epsc = exp.(range(log(150.0), log(eps[end]-100.0), length=100))
M = zeros(length(epsc))
R = zeros(length(epsc))
for i in 1:length(epsc)
    out = getMR(eps, p, n, epsc[i])
    M[i] = out[1]
    R[i] = out[2] 
end

df = DataFrame(epsc=epsc, M=M, R=R)
CSV.write("MR.csv", df, writeheader=false)