module MuTuner

using DocStringExtensions

export MuTunerLogger, update!, replay
export update_forgetful_mean
export update_forgetful_var
export update_forgetful_mv

"""
    MuTunerLogger{T<:AbstractFloat}

Type to assist with dynamically tuning the chemical potential
``\\mu`` to achieve a target particle density ``\\langle n \\rangle``.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MuTunerLogger{T<:AbstractFloat}

    "Target filling as intensive value ``\\langle n \\rangle``."
    n₀::T

    "Target filling as exensive quantity ``\\langle N \\rangle = V \\times \\langle n \\rangle``."
    N₀::T

    "Inverse temperature."
    β::T

    "System size."
    V::Int

    "Intensive energy scale."
    U₀::T

    "Extensive energy scale ``\\alpha = N / U_0``."
    α::T

    "Fraction of history to discard when calculating averages."
    c::T

    "The next chemical potential value for time `t+1`."
    μₜ₊₁::T

    "The forgetful mean of the chemical potential
    ``\\bar{\\mu}_t = \\frac{1}{L_t} \\sum_{t' = \\lceil t/c \\rceil}^t \\mu_{t'}``
    where ``L_t = t - \\lceil t/c \\rceil + 1``."
    μ_bar::T

    "The forgetful variance of the chemical potential."
    μ_var::T

    "The forgetful mean of the total particle number
    ``\\bar{N}_t = \\frac{1}{L_t} \\sum_{t' = \\lceil t/c \\rceil}^t N_{t'}``."
    N_bar::T

    "The forgetful variance of the total particle number."
    N_var::T

    "The forgetful mean of the square of total particle number
    ``\\overline{N^2}_t = \\frac{1}{L_t} \\sum_{t' = \\lceil t/c \\rceil}^t N^2_{t'}``."
    N²_bar::T

    "The forgetful mean of the compressibility
    ``\\bar{\\kappa}_t = \\frac{1}{L_t} \\sum_{t' = \\lceil t/c \\rceil}^t \\kappa_{t'}``."
    κ_bar::T

    "Sequence ``\\mu_t`` of chemical potential values."
    μ_traj::Vector{T}

    "Sequence ``N_t`` of total particle number values."
    N_traj::Vector{T}
    
    "Sequence ``N^2_t`` of total particle number squared values."
    N²_traj::Vector{T}
end

"""
    MuTunerLogger(n₀::T, β::T, V::Int, U₀::T=1.0, c::T=0.5, μ₀::T=0.0, record::Bool=false)
        where {T<:AbstractFloat}

Constructs an instance of [`MuTunerLogger`](@ref).

# Arguments
- `n₀::T`: Target particle density ``\\langle n \\rangle``.
- `β::T`: Inverse temperature.
- `V::Int`: System size.
- `U₀::T=1.0`: Characteristic intensive energy scale.
- `c::T=0.5`: Fraction of history to discard when calculating forgetful averages and variances.
- `μ₀::T=0.0`: Initial guess for chemical potential.
"""
function MuTunerLogger(n₀::T, β::T, V::Int, U₀::T=1.0, μ₀::T=0.0, c::T=0.5) where {T<:AbstractFloat}

    @assert 0.0 <= c < 1.0
    @assert V > 0

    N₀      = V*n₀
    α       = V/U₀
    μₜ₊₁    = μ₀
    μ_bar   = μ₀
    μ_var   = 0.0
    N_bar   = 0.0
    N_var   = 0.0
    N²_bar  = 0.0
    κ_bar   = α
    μ_traj  = T[]
    N_traj  = T[]
    N²_traj = T[]

    return MuTunerLogger(n₀, N₀, β, V, U₀, α, c, μₜ₊₁, μ_bar, μ_var, N_bar, N_var, N²_bar, κ_bar, μ_traj, N_traj, N²_traj)
end

MuTunerLogger(; n₀, β, V, U₀=1.0, μ₀=0, c=0.5) = MuTunerLogger(n₀, β, V, U₀, μ₀, c)


"""
    replay(μtuner::MuTunerLogger{T}) where {T<:AbstractFloat}

Replay the chemical potential tuner to its current state, returning the time series
for all relevant quantities. Note that this functions allocates all new arrays when
it returns the time series for various quantities.
"""
function replay(μtuner::MuTunerLogger{T}) where {T<:AbstractFloat}
    
    # copy existing trajectories
    μ_traj  = copy(μtuner.μ_traj)
    N_traj  = copy(μtuner.N_traj)
    N²_traj = copy(μtuner.N²_traj)

    # get current time step
    t_max = length(μ_traj) - 1

    # allocate arrays for all quantities to return
    μ_bar_traj  = zeros(T,t_max+1)
    μ_var_traj  = zeros(T,t_max+1)
    N_bar_traj  = zeros(T,t_max+1)
    N_var_traj  = zeros(T,t_max+1)
    N²_bar_traj = zeros(T,t_max+1)
    κ_bar_traj  = zeros(T,t_max+1)
    κ_fluc_traj = zeros(T,t_max+1)
    κ_min_traj  = zeros(T,t_max+1)
    κ_max_traj  = zeros(T,t_max+1)

    # iterate over time steps
    for t in 0:t_max

        # update the chemical potential based on the latest measurements
        (μₜ₊₁, μtuner.μ_bar, μtuner.μ_var, μtuner.N_bar, μtuner.N_var, μtuner.N²_bar,
        μtuner.κ_bar, κ_fluc, κ_min, κ_max) = _update!(μtuner, t)

        # record values
        μ_bar_traj[t+1]  = μtuner.μ_bar
        μ_var_traj[t+1]  = μtuner.μ_var
        N_bar_traj[t+1]  = μtuner.N_bar
        N_var_traj[t+1]  = μtuner.N_var
        N²_bar_traj[t+1] = μtuner.N²_bar
        κ_bar_traj[t+1]  = μtuner.κ_bar
        κ_fluc_traj[t+1] = κ_fluc
        κ_min_traj[t+1]  = κ_min
        κ_max_traj[t+1]  = κ_max
    end

    return (μ_traj      = μ_traj,      N_traj      = N_traj,      N²_traj    = N²_traj,
            μ_bar_traj  = μ_bar_traj,  μ_var_traj  = μ_var_traj,  N_bar_traj = N_bar_traj,
            N_var_traj  = N_var_traj,  N²_bar_traj = N²_bar_traj, κ_bar_traj = κ_bar_traj,
            κ_fluc_traj = κ_fluc_traj, κ_min_traj  = κ_min_traj,  κ_max_traj = κ_max_traj)
end


"""
    update!(μtuner::MuTunerLogger{T}, n::T, N²::T) where {T<:AbstractFloat}

Update the chemical potential given new measurements of the particle density `n`,
and the total particle number squared `N²`.
"""
function update!(μtuner::MuTunerLogger{T}, n::T, N²::T) where {T<:AbstractFloat}

    (; μₜ₊₁, μ_traj, N_traj, N²_traj, V) = μtuner

    # calculate ⟨N⟩ = V × ⟨n⟩
    N = V * n

    # record latest values for μ, ⟨N⟩ and ⟨N²⟩
    push!(μ_traj, μₜ₊₁)
    push!(N_traj, N)
    push!(N²_traj, N²)

    # get time step, keeping in mind that it is assumed to index from 0
    t = length(μ_traj) - 1

    # update the chemical potential based on the latest measurements
    μₜ₊₁, μ_bar, μ_var, N_bar, N_var, N²_bar, κ_bar, κ_fluc, κ_min, κ_max = _update!(μtuner, t)

    # record computed values as needed
    μtuner.μₜ₊₁   = μₜ₊₁
    μtuner.μ_bar  = μ_bar
    μtuner.μ_var  = μ_var
    μtuner.N_bar  = N_bar
    μtuner.N_var  = N_var
    μtuner.N²_bar = N²_bar
    μtuner.κ_bar  = κ_bar

    return μₜ₊₁
end

update!(; μtuner, n, N²) = update!(μtuner, n, N²)

# function where update to chemical potential is actually computed,
# but this method is not an exported function
function _update!(μtuner::MuTunerLogger{T}, t::Int) where {T<:AbstractFloat}

    (; μ_bar, μ_var, κ_bar, N_bar, N_var, N²_bar, N₀, β, α, c) = μtuner

    # update forgetful mean and variance of μ
    μ_traj = @view μtuner.μ_traj[1:t+1]
    μ_bar, μ_var = update_forgetful_mv(μ_traj, μ_bar, μ_var, c)

    # update forgetful mean of ⟨N⟩
    N_traj = @view μtuner.N_traj[1:t+1]
    N_bar = update_forgetful_mean(N_traj, N_bar, c)

    # update forgetful mean of ⟨N²⟩
    N²_traj = @view μtuner.N²_traj[1:t+1]
    N²_bar = update_forgetful_mean(N²_traj, N²_bar, c)

    # update forgetful variance of ⟨N⟩
    N_var = N²_bar - N_bar^2

    # get fluctuation based estimate of κ
    κ_fluc = β * N_var

    # get lower bound for κ
    κ_min = α/sqrt(t+1)

    # get upper bound for κ
    if iszero(t) || iszero(μ_var)
        κ_max = κ_min
    else
        κ_max = sqrt(N_var / μ_var)
    end

    # calculate bounded estimate of κ
    κ_bar = max(κ_min , min(κ_max , κ_fluc))

    # update the chemical potential
    μₜ₊₁ = μ_bar + (N₀ - N_bar)/κ_bar

    return (μₜ₊₁, μ_bar, μ_var, N_bar, N_var, N²_bar, κ_bar, κ_fluc, κ_min, κ_max)
end


"""
    update_forgetful_mean(x::AbstractVector{T}, x̄ₜ::T, c::T)
        where {T<:AbstractFloat}

Given the previous value of the forgetful mean `x̄ₜ`, calculate its updated value
`x̄ₜ₊₁` assuming that `x[end] = xₜ₊₁` has already been appended to `x`.
The oldest `c` fraction of values is discarded when calculating the forgetful mean.
"""
function update_forgetful_mean(x::AbstractVector{T}, x̄ₜ::T, c::T) where {T<:AbstractFloat}

    @assert 0.0 <= c < 1.0

    # t is assumed to index from 0 for algorithm, but arrays index from 1.
    # keep in mind when reading the code below.

    if length(x) == 1
        return x[1]
    end

    # assume xₜ₊₁ = x[end].
    # therefore, (t + 2) = length(x) ==> (t + 1) = (length(x) - 1).
    tp1  = length(x) - 1
    t    = tp1 - 1
    ct   = ceil(Int,c*t)
    ctp1 = ceil(Int,c*tp1)
    Lₜ₊₁ = tp1 - ctp1 + 1
    if ctp1 == ct
        x̄ₜ₊₁ = x̄ₜ + (x[tp1+1] - x̄ₜ)/Lₜ₊₁
    else
        x̄ₜ₊₁ = x̄ₜ + (x[tp1+1] - x[ct+1])/Lₜ₊₁
    end

    return x̄ₜ₊₁
end

update_forgetful_mean(; x, x̄ₜ, c) = update_forgetful_mean(x, x̄ₜ, c)


"""
    update_forgetful_var(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ₊₁::T, x̄ₜ::T, c::T)
        where {T<:AbstractFloat}

Given the previous value for the forgetful variance `V̄ₜ` and foregetul mean `x̄ₜ`,
and the updated value for the forgetful mean `x̄ₜ₊₁`, calculate the updated
value `V̄ₜ₊₁`, assuming that `x[end] = xₜ₊₁` has already been appended to `x`.
The oldest `c` fraction of values is discarded when calculating the forgetful variance.
"""
function update_forgetful_var(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ₊₁::T, x̄ₜ::T, c::T) where {T<:AbstractFloat}

    @assert 0.0 <= c < 1.0

    # t is assumed to index from 0 for algorithm, but arrays index from 1.
    # keep in mind when reading the code below.

    if length(x) == 1
        return 0.0
    end

    # assume xₜ₊₁ = x[end].
    # therefore, (t + 2) = length(x) ==> (t + 1) = (length(x) - 1).
    tp1  = length(x) - 1
    t    = tp1 - 1
    ct   = ceil(Int,c*t)
    ctp1 = ceil(Int,c*tp1)
    Lₜ   = t - ct + 1
    Mₜ   = V̄ₜ * Lₜ
    if ctp1 == ct
        Mₜ₊₁ = Mₜ + (x[tp1+1]-x̄ₜ)*(x[tp1+1]-x̄ₜ₊₁)
    else
        Mₜ₊₁ = Mₜ + (x[tp1+1]-x[ct+1])*(x[tp1+1]-x̄ₜ₊₁+x[ct+1]-x̄ₜ)
    end
    Lₜ₊₁ = tp1 - ctp1 + 1
    V̄ₜ₊₁ = Mₜ₊₁/Lₜ₊₁

    return V̄ₜ₊₁
end

update_forgetful_var(; x, V̄ₜ, x̄ₜ₊₁, x̄ₜ, c) = update_forgetful_var(x, V̄ₜ, x̄ₜ₊₁, x̄ₜ, c)


"""
    update_forgetful_mv(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ::T, c::T)
        where {T<:AbstractFloat}

Given the previous value of the forgetful variance `V̄ₜ` and foregetul mean `x̄ₜ`, calculate their updated
values `V̄ₜ₊₁` and `x̄ₜ₊₁`, assuming that `x[end] = xₜ₊₁` has already been appended to `x`.
The oldest `c` fraction of values is discarded when calculating the forgetful mean and variance.
"""
function update_forgetful_mv(x::AbstractVector{T}, x̄ₜ::T, V̄ₜ::T, c::T) where {T<:AbstractFloat}

    @assert 0.0 <= c < 1.0

    x̄ₜ₊₁ = update_forgetful_mean(x, x̄ₜ, c)
    V̄ₜ₊₁ = update_forgetful_var(x, V̄ₜ, x̄ₜ₊₁, x̄ₜ, c)

    return x̄ₜ₊₁, V̄ₜ₊₁
end

update_forgetful_mv(; x, x̄ₜ, V̄ₜ, c) = update_forgetful_mv(x, x̄ₜ, V̄ₜ, c)

end