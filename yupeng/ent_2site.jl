using DelimitedFiles, LinearAlgebra, LindbladEq
#--------------------------------------------------------------------------------
"""
H = ∑ cᵢ⁺cᵢ₊₁ + cᵢ₊₁⁺cᵢ.
"""
function hopping_ham(L::Integer; PBC::Bool=false)
    H = diagm(1 => ones(L-1), -1 => ones(L-1))
    PBC && (H[L,1] = H[1,L] = 1)
    H
end
#--------------------------------------------------------------------------------
"""
Jump operator:
    Lᵢ = √γ * dᵢ⁺dᵢ,
    dᵢ⁺ = (cᵢ₋₁⁺ + 1im cᵢ⁺ ) / √2.
"""
function nodal_jump(L::Integer;)
    V = [1 , 1im] / sqrt(2)
    l1 = [QuasiMode([i, mod(i,L)+1], V, L) for i in 1:2:L-1]
    l2 = [QuasiMode([i, mod(i,L)+1], V, L) for i in 2:2:L]
    l1, l2
end

#--------------------------------------------------------------------------------
# Main
#--------------------------------------------------------------------------------
function ent_evo(;L=128, γ=1.0, dt=0.05, T=500, a=1, config="Z2", p::Integer=20)
    N = round(Int, T/dt)
    evo = evo_operator(Hermitian(hopping_ham(L, PBC=true)), dt)
    l1, l2 = nodal_jump(L; a)
    s = FreeFermionState(L=L, N=L÷2, config=config)
    EE = zeros(N÷p)
    hL = L÷2
    for k = 1:N
        s = evo * s
        wiener!(l1, s, γ*dt)
        wiener!(l2, s, γ*dt)
        if iszero(mod(k,p))
            ee = zeros(hL)
            Threads.@threads for i in 1:hL 
                ind = mod.(i-1:hL+i-2, L) .+ 1
                ee[i] = ent_S(s, ind)
            end
            EE[k÷p] = sum(ee) / hL
            writedlm("$L.dat", EE)
            GC.gc()
        end
    end
end

Ls = round.(Int, 10 .^ (1:0.1:2.5)) .÷ 2 .* 2
for L in Ls
    ent_evo(;L)
end

