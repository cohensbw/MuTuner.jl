module MuTuner

using Printf

export MuTunerLogger, update!, save, replay
export update_forgetful_mean
export update_forgetful_var
export update_forgetful_mv

@doc raw"""
    MuTunerLogger{T<:AbstractFloat, S<:Number}

Type to assist with dynamically tuning the chemical potential
``mu`` to achieve a target particle density ``\langle n \rangle = n_0``.

# Fields

- `n₀::T`: Target filling as intensive value ``\langle n \rangle = n_0``.
- `N₀::T`: Target filling as exensive quantity ``\langle N \rangle = N_0 = V \cdot n_0``.
- `β::T`: Inverse temperature.
- `V::Int`: System size.
- `u₀::T`: Intensive energy scale.
- `α::T`: Extensive energy scale ``\alpha = N / U_0``.
- `c::T`: Fraction of history to discard when calculating averages.
- `μ_tp1::T`: Next chemical potential value for time `t+1`.
- `μ_bar::T`: Forgetful mean of the chemical potential.
- `μ_var::T`: Forgetful variance of the chemical potential.
- `N_bar::S`: Forgetful mean of the total particle number, accounting for the sign.
- `N_var::T`: Forgetful variance of the total particle number, accounting for the sign.
- `s_bar::T`: Forgetful mean of the sign.
- `s_var::T`: Forgetful variance of the sign.
- `N²_bar::S`: Forgetful mean of the square of total particle number, accounting for the sign.
- `κ_bar::T`: Forgetful mean of the compressibility.
- `μ_traj::Vector{T}`: Timeseries of chemical potential values, ``\mu_t.``
- `N_traj::Vector{S}`: Timeseries of total particle number values accounting for the sign, ``s_t \cdot N_t.``
- `s_traj::Vector{S}`: Timeseries of the sign, ``s_t.``
- `N²_traj::Vector{S}`: Timeseries of total particle number square values accounting for the sign, ``s_t \cdot N^2_t.``
"""
mutable struct MuTunerLogger{T<:AbstractFloat, S<:Number}

    "Target filling as intensive value."
    n₀::T

    "Target filling as exensive quantity."
    N₀::T

    "Inverse temperature."
    β::T

    "System size."
    V::Int

    "Intensive energy scale."
    u₀::T

    "Extensive energy scale."
    α::T

    "Fraction of history to discard when calculating averages."
    c::T

    "The next chemical potential value for time `t+1`."
    μ_tp1::T

    "The forgetful mean of the chemical potential"
    μ_bar::T

    "The forgetful variance of the chemical potential."
    μ_var::T

    "The forgetful mean of the total particle number"
    N_bar::S

    "The forgetful variance of the total particle number."
    N_var::T

    "The forgetful mean of the sign."
    s_bar::S

    "The forgetful variance of the sign."
    s_var::S

    "The forgetful mean of the square of total particle number."
    N²_bar::S

    "The forgetful mean of the compressibility."
    κ_bar::T

    "Sequence of chemical potential values."
    μ_traj::Vector{T}

    "Sequence of total particle number values."
    N_traj::Vector{S}

    "Sequence of signs."
    s_traj::Vector{S}
    
    "Sequence of total particle number squared values."
    N²_traj::Vector{S}
end

@doc raw"""
    MuTunerLogger(n₀::T, β::T, V::Int, u₀::T=1.0, μ₀::T=0.0, c::T=0.5,
                  s::S=zero(T)) where {T<:AbstractFloat, S<:Number}

Constructs an instance of [`MuTunerLogger`](@ref).

# Arguments
- `n₀::T`: Target particle density ``\langle n \rangle``.
- `β::T`: Inverse temperature.
- `V::Int`: System size.
- `u₀::T=1.0`: Characteristic intensive energy scale.
- `μ₀::T=0.0`: Initial guess for chemical potential.
- `c::T=0.5`: Fraction of history to discard when calculating forgetful averages and variances.
- `s::S=zero(T)`: An example of the sign, to determine whether the type of the sign if real or complex.
"""
function MuTunerLogger(n₀::T, β::T, V::Int, u₀::T=1.0, μ₀::T=0.0, c::T=0.5, s::S=zero(T)) where {T<:AbstractFloat, S<:Number}

    @assert 0.0 <= c < 1.0
    @assert V > 0

    N₀      = V*n₀
    α       = V/u₀
    μ_tp1   = μ₀
    μ_bar   = μ₀
    μ_var   = zero(T)
    N_bar   = zero(S)
    N_var   = zero(T)
    s_bar   = zero(S)
    s_var   = zero(S)
    N²_bar  = zero(S)
    κ_bar   = α
    μ_traj  = T[]
    N_traj  = S[]
    s_traj  = S[]
    N²_traj = S[]

    return MuTunerLogger(n₀, N₀, β, V, u₀, α, c, μ_tp1, μ_bar, μ_var, N_bar, N_var, s_bar, s_var, N²_bar, κ_bar, μ_traj, N_traj, s_traj, N²_traj)
end

MuTunerLogger(; n₀, β, V, u₀=1.0, μ₀=0.0, c=0.5, s=1.0) = MuTunerLogger(n₀, β, V, u₀, μ₀, c, s)


@doc raw"""
    save(μtuner::MuTunerLogger{T,S}, filename::String, filepath::String = "")

Replay the chemical potential tuner to its current state, writing the time series of
relevant values to a space-delimited CSV file.
"""
function save(μtuner::MuTunerLogger{T,S}, filename::String, filepath::String = "") where {T<:AbstractFloat, S<:Number}

    # write csv file
    open(joinpath(filepath,filename), "w") do fout
        # write file header
        write(fout, "t mu N Nsqrd sign_real sign_imag mu_bar mu_var N_bar N_var Nsqrd_bar kappa_bar kappa_fluc kappa_min kappa_max sign_bar_real sign_bar_imag sign_var_real sign_var_imag\n")
        # get current of time-steps
        t_max = length(μtuner.μ_traj) - 1
        # iterate over time steps
        for t in 0:t_max
            # update the chemical potential based on the latest measurements
            (μ_tp1, μtuner.μ_bar, μtuner.μ_var, μtuner.N_bar, μtuner.N_var, μtuner.s_bar, μtuner.s_var,
            μtuner.N²_bar, μtuner.κ_bar, κ_fluc, κ_min, κ_max) = _update!(μtuner, t)
            # write data to file
            @printf(fout, "%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
                    t, real(μtuner.μ_traj[t+1]), real(μtuner.N_traj[t+1]), real(μtuner.N²_traj[t+1]), real(μtuner.s_traj[t+1]), imag(μtuner.s_traj[t+1]),
                    real(μtuner.μ_bar), real(μtuner.μ_var), real(μtuner.N_bar), real(μtuner.N_var), real(μtuner.N²_bar), real(μtuner.κ_bar),
                    real(κ_fluc), real(κ_min), real(κ_max), real(μtuner.s_bar), imag(μtuner.s_bar), real(μtuner.s_var), imag(μtuner.s_var))
        end
    end

    return nothing
end


@doc raw"""
    replay(μtuner::MuTunerLogger{T}) where {T<:AbstractFloat}

Replay the chemical potential tuner to its current state, returning the time series
for all relevant quantities. Note that this functions allocates all new arrays when
it returns the time series for various quantities.
"""
function replay(μtuner::MuTunerLogger{T,S}) where {T<:AbstractFloat, S<:Number}
    
    # copy existing trajectories
    μ_traj  = copy(μtuner.μ_traj)
    N_traj  = copy(μtuner.N_traj)
    s_traj  = copy(μtuner.s_traj)
    N²_traj = copy(μtuner.N²_traj)

    # get current time step
    t_max = length(μ_traj) - 1

    # allocate arrays for all quantities to return
    μ_bar_traj  = zeros(T,t_max+1)
    μ_var_traj  = zeros(T,t_max+1)
    N_bar_traj  = zeros(T,t_max+1)
    N_var_traj  = zeros(T,t_max+1)
    s_bar_traj  = zeros(S,t_max+1)
    s_var_traj  = zeros(S,t_max+1)
    N²_bar_traj = zeros(T,t_max+1)
    κ_bar_traj  = zeros(T,t_max+1)
    κ_fluc_traj = zeros(T,t_max+1)
    κ_min_traj  = zeros(T,t_max+1)
    κ_max_traj  = zeros(T,t_max+1)

    # iterate over time steps
    for t in 0:t_max

        # update the chemical potential based on the latest measurements
        (μ_tp1, μtuner.μ_bar, μtuner.μ_var, μtuner.N_bar, μtuner.N_var, μtuner.s_bar, μtuner.s_var,
        μtuner.N²_bar, μtuner.κ_bar, κ_fluc, κ_min, κ_max) = _update!(μtuner, t)

        # record values
        μ_bar_traj[t+1]  = μtuner.μ_bar
        μ_var_traj[t+1]  = μtuner.μ_var
        N_bar_traj[t+1]  = μtuner.N_bar
        N_var_traj[t+1]  = μtuner.N_var
        s_bar_traj[t+1]  = μtuner.s_bar
        s_var_traj[t+1]  = μtuner.s_var
        N²_bar_traj[t+1] = μtuner.N²_bar
        κ_bar_traj[t+1]  = μtuner.κ_bar
        κ_fluc_traj[t+1] = κ_fluc
        κ_min_traj[t+1]  = κ_min
        κ_max_traj[t+1]  = κ_max
    end

    return (μ_traj      = μ_traj,      N_traj      = N_traj,      N²_traj    = N²_traj,
            μ_bar_traj  = μ_bar_traj,  μ_var_traj  = μ_var_traj,  N_bar_traj = N_bar_traj,
            N_var_traj  = N_var_traj,  N²_bar_traj = N²_bar_traj, κ_bar_traj = κ_bar_traj,
            κ_fluc_traj = κ_fluc_traj, κ_min_traj  = κ_min_traj,  κ_max_traj = κ_max_traj,
            s_bar_traj  = s_bar_traj,  s_var_traj  = s_var_traj,  s_traj     = s_traj)
end


@doc raw"""
    update!(μtuner::MuTunerLogger{T}, n::T, N²::T,
            s::S=one(S)) where {T<:AbstractFloat, S<:Number}

Update the chemical potential given new measurements of the particle density `n`,
the total particle number squared `N²`, and the sign `s`.
"""
function update!(μtuner::MuTunerLogger{T}, n::T, N²::T, s::S=one(S)) where {T<:AbstractFloat, S<:Number}

    (; μ_tp1, μ_traj, N_traj, s_traj, N²_traj, V) = μtuner

    # calculate ⟨N⟩ = V × ⟨n⟩
    N = V * n

    # record latest values for μ, ⟨s⋅N⟩, s and ⟨s⋅N²⟩
    push!(μ_traj, μ_tp1)
    push!(s_traj, s)
    push!(N_traj, s*N)
    push!(N²_traj, s*N²)

    # get time step, keeping in mind that it is assumed to index from 0
    t = length(μ_traj) - 1

    # update the chemical potential based on the latest measurements
    μ_tp1, μ_bar, μ_var, N_bar, N_var, s_bar, s_var, N²_bar, κ_bar, κ_fluc, κ_min, κ_max = _update!(μtuner, t)

    # record computed values as needed
    μtuner.μ_tp1   = μ_tp1
    μtuner.μ_bar  = μ_bar
    μtuner.μ_var  = μ_var
    μtuner.N_bar  = N_bar
    μtuner.N_var  = N_var
    μtuner.s_bar  = s_bar
    μtuner.s_var  = s_var
    μtuner.N²_bar = N²_bar
    μtuner.κ_bar  = κ_bar

    return μ_tp1
end

update!(; μtuner, n, N², s=1.0) = update!(μtuner, n, N², s)

# function where update to chemical potential is actually computed,
# but this method is not an exported function
function _update!(μtuner::MuTunerLogger{T,S}, t::Int) where {T<:AbstractFloat, S<:Number}

    (; μ_bar, μ_var, κ_bar, N_bar, N_var, s_bar, s_var, N²_bar, N₀, β, α, c) = μtuner

    # update forgetful mean and variance of μ
    μ_traj = @view μtuner.μ_traj[1:t+1]
    μ_bar, μ_var = update_forgetful_mv(μ_traj, μ_bar, μ_var, c)

    # update forgetful mean of ⟨s⋅N⟩
    N_traj = @view μtuner.N_traj[1:t+1]
    N_bar = update_forgetful_mean(N_traj, N_bar, c)

    # update forgetful mean of ⟨s⋅N²⟩
    N²_traj = @view μtuner.N²_traj[1:t+1]
    N²_bar = update_forgetful_mean(N²_traj, N²_bar, c)

    # update forgetful mean of the sign ⟨s⟩
    s_traj = @view μtuner.s_traj[1:t+1]
    s_bar, s_var = update_forgetful_mv(s_traj, s_bar, s_var, c)

    # update forgetful variance of based on ⟨N²⟩-⟨N⟩², accounting for
    # the average sign and reweighting ⟨O⟩=⟨s⋅O⟩/⟨s⟩
    if abs(s_bar) ≈ 0
        N_var = T(Inf)
    else
        N_var = real(N²_bar/s_bar - (N_bar/s_bar)^2)
    end

    # get fluctuation based estimate of κ
    κ_fluc = β * N_var

    # get lower bound for κ
    κ_min = α/sqrt(t+1)

    # get upper bound for κ
    if iszero(t) || iszero(μ_var) || (N_var / μ_var) < 0
        κ_max = κ_min
    else
        κ_max = sqrt(N_var / μ_var)
    end

    # calculate bounded estimate of κ
    κ_bar = max(κ_min , min(κ_max , κ_fluc))

    # update the chemical potential
    μ_tp1 = μ_bar + real(N₀ - (N_bar/s_bar))/κ_bar

    return (μ_tp1, μ_bar, μ_var, N_bar, N_var, s_bar, s_var, N²_bar, κ_bar, κ_fluc, κ_min, κ_max)
end


@doc raw"""
    update_forgetful_mean(x::AbstractVector{T}, x̄ₜ::T, c::E) where {T<:Number, E<:AbstractFloat}

Given the previous value of the forgetful mean `x̄ₜ`, calculate its updated value
`x̄ₜ₊₁` assuming that `x[end] = xₜ₊₁` has already been appended to `x`.
The oldest `c` fraction of values is discarded when calculating the forgetful mean.
"""
function update_forgetful_mean(x::AbstractVector{T}, x̄ₜ::T, c::E) where {T<:Number, E<:AbstractFloat}

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
    update_forgetful_var(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ₊₁::T, x̄ₜ::T, c::E) where {T<:Number, E<:AbstractFloat}

Given the previous value for the forgetful variance `V̄ₜ` and foregetul mean `x̄ₜ`,
and the updated value for the forgetful mean `x̄ₜ₊₁`, calculate the updated
value `V̄ₜ₊₁`, assuming that `x[end] = xₜ₊₁` has already been appended to `x`.
The oldest `c` fraction of values is discarded when calculating the forgetful variance.
"""
function update_forgetful_var(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ₊₁::T, x̄ₜ::T, c::E) where {T<:Number, E<:AbstractFloat}

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


@doc raw"""
    update_forgetful_mv(x::AbstractVector{T}, V̄ₜ::T, x̄ₜ::T, c::E) where {T<:Number, E<:AbstractFloat}

Given the previous value of the forgetful variance `V̄ₜ` and foregetul mean `x̄ₜ`, calculate their updated
values `V̄ₜ₊₁` and `x̄ₜ₊₁`, assuming that `x[end] = xₜ₊₁` has already been appended to `x`.
The oldest `c` fraction of values is discarded when calculating the forgetful mean and variance.
"""
function update_forgetful_mv(x::AbstractVector{T}, x̄ₜ::T, V̄ₜ::T, c::E) where {T<:Number, E<:AbstractFloat}

    @assert 0.0 <= c < 1.0

    x̄ₜ₊₁ = update_forgetful_mean(x, x̄ₜ, c)
    V̄ₜ₊₁ = update_forgetful_var(x, V̄ₜ, x̄ₜ₊₁, x̄ₜ, c)

    return x̄ₜ₊₁, V̄ₜ₊₁
end

update_forgetful_mv(; x, x̄ₜ, V̄ₜ, c) = update_forgetful_mv(x, x̄ₜ, V̄ₜ, c)

end