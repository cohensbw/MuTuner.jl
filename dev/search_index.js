var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [MuTuner]","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [MuTuner]","category":"page"},{"location":"api/#MuTuner.MuTunerLogger","page":"API","title":"MuTuner.MuTunerLogger","text":"MuTunerLogger{T<:AbstractFloat}\n\nType to assist with dynamically tuning the chemical potential mu to achieve a target particle density langle n rangle.\n\nFields\n\nn₀::AbstractFloat\nTarget filling as intensive value langle n rangle.\nN₀::AbstractFloat\nTarget filling as exensive quantity langle N rangle = V times langle n rangle.\nβ::AbstractFloat\nInverse temperature.\nV::Int64\nSystem size.\nU₀::AbstractFloat\nIntensive energy scale.\nα::AbstractFloat\nExtensive energy scale alpha = N  U_0.\nc::AbstractFloat\nFraction of history to discard when calculating averages.\nμₜ₊₁::AbstractFloat\nThe next chemical potential value for time t+1.\nμ_bar::AbstractFloat\nThe forgetful mean of the chemical potential     barmu_t = frac1L_t sum_t = lceil tc rceil^t mu_t     where L_t = t - lceil tc rceil + 1.\nμ_var::AbstractFloat\nThe forgetful variance of the chemical potential.\nN_bar::AbstractFloat\nThe forgetful mean of the total particle number     barN_t = frac1L_t sum_t = lceil tc rceil^t N_t.\nN_var::AbstractFloat\nThe forgetful variance of the total particle number.\nN²_bar::AbstractFloat\nThe forgetful mean of the square of total particle number     overlineN^2_t = frac1L_t sum_t = lceil tc rceil^t N^2_t.\nκ_bar::AbstractFloat\nThe forgetful mean of the compressibility     barkappa_t = frac1L_t sum_t = lceil tc rceil^t kappa_t.\nμ_traj::Vector{T} where T<:AbstractFloat\nSequence mu_t of chemical potential values.\nN_traj::Vector{T} where T<:AbstractFloat\nSequence N_t of total particle number values.\nN²_traj::Vector{T} where T<:AbstractFloat\nSequence N^2_t of total particle number squared values.\n\n\n\n\n\n","category":"type"},{"location":"api/#MuTuner.MuTunerLogger-Union{Tuple{T}, Tuple{T, T, Int64}, Tuple{T, T, Int64, T}, Tuple{T, T, Int64, T, T}, Tuple{T, T, Int64, T, T, T}} where T<:AbstractFloat","page":"API","title":"MuTuner.MuTunerLogger","text":"MuTunerLogger(n₀::T, β::T, V::Int, U₀::T=1.0, c::T=0.5, μ₀::T=0.0, record::Bool=false)\n    where {T<:AbstractFloat}\n\nConstructs an instance of MuTunerLogger.\n\nArguments\n\nn₀::T: Target particle density langle n rangle.\nβ::T: Inverse temperature.\nV::Int: System size.\nU₀::T=1.0: Characteristic intensive energy scale.\nc::T=0.5: Fraction of history to discard when calculating forgetful averages and variances.\nμ₀::T=0.0: Initial guess for chemical potential.\n\n\n\n\n\n","category":"method"},{"location":"api/#MuTuner.replay-Union{Tuple{MuTunerLogger{T}}, Tuple{T}} where T<:AbstractFloat","page":"API","title":"MuTuner.replay","text":"replay(μtuner::MuTunerLogger{T}) where {T<:AbstractFloat}\n\nReplay the chemical potential tuner to its current state, returning the time series for all relevant quantities. Note that this functions allocates all new arrays when it returns the time series for various quantities.\n\n\n\n\n\n","category":"method"},{"location":"api/#MuTuner.update!-Union{Tuple{T}, Tuple{MuTunerLogger{T}, T, T}} where T<:AbstractFloat","page":"API","title":"MuTuner.update!","text":"update!(μtuner::MuTunerLogger{T}, n::T, N²::T) where {T<:AbstractFloat}\n\nUpdate the chemical potential given new measurements of the particle density n, and the total particle number squared N².\n\n\n\n\n\n","category":"method"},{"location":"api/#MuTuner.update_forgetful_mean-Union{Tuple{T}, Tuple{AbstractVector{T}, T, T}} where T<:AbstractFloat","page":"API","title":"MuTuner.update_forgetful_mean","text":"update_forgetful_mean(x::AbstractVector{T}, x̄ₜ::T, c::T)\n    where {T<:AbstractFloat}\n\nGiven the previous value of the forgetful mean x̄ₜ, calculate its updated value x̄ₜ₊₁ assuming that x[end] = xₜ₊₁ has already been appended to x. The oldest c fraction of values is discarded when calculating the forgetful mean.\n\n\n\n\n\n","category":"method"},{"location":"api/#MuTuner.update_forgetful_mv-Union{Tuple{T}, Tuple{AbstractVector{T}, T, T, T}} where T<:AbstractFloat","page":"API","title":"MuTuner.update_forgetful_mv","text":"update_forgetful_mv(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ::T, c::T)\n    where {T<:AbstractFloat}\n\nGiven the previous value of the forgetful variance V̄ₜ and foregetul mean x̄ₜ, calculate their updated values V̄ₜ₊₁ and x̄ₜ₊₁, assuming that x[end] = xₜ₊₁ has already been appended to x. The oldest c fraction of values is discarded when calculating the forgetful mean and variance.\n\n\n\n\n\n","category":"method"},{"location":"api/#MuTuner.update_forgetful_var-Union{Tuple{T}, Tuple{AbstractVector{T}, T, T, T, T}} where T<:AbstractFloat","page":"API","title":"MuTuner.update_forgetful_var","text":"update_forgetful_var(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ₊₁::T, x̄ₜ::T, c::T)\n    where {T<:AbstractFloat}\n\nGiven the previous value for the forgetful variance V̄ₜ and foregetul mean x̄ₜ, and the updated value for the forgetful mean x̄ₜ₊₁, calculate the updated value V̄ₜ₊₁, assuming that x[end] = xₜ₊₁ has already been appended to x. The oldest c fraction of values is discarded when calculating the forgetful variance.\n\n\n\n\n\n","category":"method"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"On this page we include examples for how this package can be used.","category":"page"},{"location":"examples/#The-Ising-Model","page":"Examples","title":"The Ising Model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"While the algorithm implemented in this package, and introduced in Phys. Rev. E 105, 045311, was initially designed to tune the chemical potential mu to achieve a target particle density n_0, it can also be applied to achieve a target magnetization m_0 in a spin model by tuning the applied magnetic field B.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example we use MuTuner.jl to tune the magnetic field B to achieve a target magnetization m_0 in a Monte Carlo simulation of the square lattice Ising model","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"H = -J sum_langle ijrangle s_is_j - B sum_i s_i","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In that above J is the coupling strength, and each Ising field can take on a value s_i = pm 1. The tuning algorithm can be directly applied by applying the substitutions","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"mu mapsto B quad N mapsto M quad kappa mapsto chi","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"where kappa and chi are the compressibility and magnetic susceptibility respectively, and M is the total magnetization. Therefore, given a fixed temperature T, we are interested in tuning the average magnetization langle m rangle = langle M rangle  L^2 to the target magnetization m_0, where L is the linear extent of the square lattice.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Below is an example script using MuTuner.jl in a Monte Carlo simulation of the Ising Model:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using MuTuner\n\n# calculate the site energy for site (i,j)\nfunction calc_Eᵢⱼ(s::Matrix{Int}, i::Int, j::Int, J::E, B::E) where {E<:AbstractFloat}\n    \n    L    = size(s,1)\n    Eᵢⱼ  = -B * s[i,j]\n    Eᵢⱼ += -J * s[i,j] * s[mod1(i+1,L),j]\n    Eᵢⱼ += -J * s[i,j] * s[mod1(i-1,L),j]\n    Eᵢⱼ += -J * s[i,j] * s[i,mod1(j+1,L)]\n    Eᵢⱼ += -J * s[i,j] * s[i,mod1(j-1,L)]\n    return Eᵢⱼ\nend\n\n# calculate the change in energy associated with a spin flip on site (i,j)\nfunction calc_ΔEᵢⱼ(s::Matrix{Int}, i::Int, j::Int, J::E, B::E) where {E<:AbstractFloat}\n    \n    Eᵢⱼ  = calc_Eᵢⱼ(s,i,j,J,B)\n    ΔEᵢⱼ = -2*Eᵢⱼ\n    return ΔEᵢⱼ\nend\n\n# sweep through the lattice, attempting a spin flip on each site in the lattice\nfunction monte_carlo_sweep!(s::Matrix{Int}, T::E, J::E, B::E) where {E<:AbstractFloat}\n    \n    # get lattice size\n    L = size(s,1)\n    \n    # iterate over all sites in lattice\n    for j in 1:L\n        for i in 1:L\n            \n            # calculate change in energy associated with spin flip\n            ΔEᵢⱼ = calc_ΔEᵢⱼ(s,i,j,J,B)\n            \n            # calculate metropolis-hastings acceptance probability\n            Pᵢⱼ  = min(1.0, exp(-ΔEᵢⱼ/T))\n            \n            # accept or reject proposed spin flip\n            if rand() < Pᵢⱼ\n                s[i,j] = -s[i,j]\n            end\n        end\n    end\n    \n    return nothing\nend\n\n# system size\nL = 100\n\n# temperature\nT = 2.5\n\n# coupling strength\nJ = 1.0\n\n# initial magnetic field\nB = 0.0\n\n# number of Monte Carlo sweeps to perform\nN = 100_000\n\n# target magnetization\nm₀ = 0.5\n\n# initialize ising spins to random initial configuration\ns = rand(-1:2:1,L,L)\n\n# initialize MuTunerLogger for tuning magnetic field B\nBtuner = MuTunerLogger(n₀=m₀, β=1/T, V=L^2, U₀=J, μ₀=B, c=0.5)\n\n# perform Monte Carlo simulation\nfor i in 1:N\n    \n    # perform Monte Carlo sweep through lattice\n    monte_carlo_sweep!(s, T, J, B)\n    \n    # measure ⟨m⟩ = ⟨M⟩/L²\n    M = float(sum(s))\n    m = M/L^2\n    \n    # measure ⟨M²⟩\n    M² = M^2\n    \n    # update the magnetic field\n    B = update!(μtuner=Btuner, n=m, N²=M²)\nend\n\n# replay the tuning time-series for all relevant observables.\n# these resulting time-series arrays can be used to generate a figure\n# comparable to Figure 1 in the paper Phys. Rev. E 105, 045311.\ntape        = replay(Btuner);\nB_traj      = tape.μ_traj\nB_bar_traj  = tape.μ_bar_traj\nm_traj      = tape.N_traj / L^2\nm_bar_traj  = tape.N_bar_traj / L^2\nχ_bar_traj  = tape.κ_bar_traj\nχ_fluc_traj = tape.κ_fluc_traj\nχ_min_traj  = tape.κ_min_traj\nχ_max_traj  = tape.κ_max_traj","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MuTuner","category":"page"},{"location":"#MuTuner.jl","page":"Home","title":"MuTuner.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MuTuner.jl. This package exports a method for tuning the chemical potential to achieve a target particle number in grand canonical Monte Carlo simulations. The algorithm implemented in this package was introduced in the paper Phys. Rev. E 105, 045311. If you use this package in work leading to a publication, please consider citing this paper.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install MuTuner.jl run following in the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add https://github.com/cohensbw/MuTuner.jl","category":"page"}]
}
