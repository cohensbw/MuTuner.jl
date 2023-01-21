# Examples

On this page we include examples for how this package can be used.

## The Ising Model

While the algorithm implemented in this package, and introduced in
[Phys. Rev. E 105, 045311](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.105.045311),
was initially designed to tune the chemical potential ``\mu`` to achieve a target particle
density ``n_0``, it can also be applied to achieve a target magnetization
``m_0`` in a spin model by tuning the applied magnetic field ``B``.

In this example we use [`MuTuner.jl`](https://github.com/cohensbw/MuTuner.jl) to tune the
magnetic field ``B`` to achieve a target magnetization ``m_0`` in a Monte Carlo
simulation of the square lattice Ising model
```math
H = -J \sum_{\langle ij\rangle} s_{i}s_{j} - B \sum_{i} s_{i},
```
In that above ``J`` is the coupling strength, and each
Ising field can take on a value ``s_i = \pm 1``.
The tuning algorithm can be directly applied by applying the substitutions
```math
\mu \mapsto B, \quad N \mapsto M, \quad \kappa \mapsto \chi,
```
where ``\kappa`` and ``\chi`` are the compressibility and magnetic
susceptibility respectively, and ``M`` is the total magnetization.
Therefore, given a fixed temperature ``T``, we are interested in tuning the average
magnetization ``\langle m \rangle = \langle M \rangle / L^2`` to the target magnetization
``m_0``, where ``L`` is the linear extent of the square lattice.


Below is an example script using [`MuTuner.jl`](https://github.com/cohensbw/MuTuner.jl) in a
Monte Carlo simulation of the Ising Model:

```julia
using Revise
using MuTuner

# calculate the site energy for site (i,j)
function calc_Eᵢⱼ(s::Matrix{Int}, i::Int, j::Int, J::E, B::E) where {E<:AbstractFloat}
    
    L    = size(s,1)
    Eᵢⱼ  = -B * s[i,j]
    Eᵢⱼ += -J * s[i,j] * s[mod1(i+1,L),j]
    Eᵢⱼ += -J * s[i,j] * s[mod1(i-1,L),j]
    Eᵢⱼ += -J * s[i,j] * s[i,mod1(j+1,L)]
    Eᵢⱼ += -J * s[i,j] * s[i,mod1(j-1,L)]
    return Eᵢⱼ
end

# calculate the change in energy associated with a spin flip on site (i,j)
function calc_ΔEᵢⱼ(s::Matrix{Int}, i::Int, j::Int, J::E, B::E) where {E<:AbstractFloat}
    
    Eᵢⱼ  = calc_Eᵢⱼ(s,i,j,J,B)
    ΔEᵢⱼ = -2*Eᵢⱼ
    return ΔEᵢⱼ
end

# sweep through the lattice, attempting a spin flip on each site in the lattice
function monte_carlo_sweep!(s::Matrix{Int}, T::E, J::E, B::E) where {E<:AbstractFloat}
    
    # get lattice size
    L = size(s,1)
    
    # iterate over all sites in lattice
    @fastmath @inbounds for j in 1:L
        for i in 1:L
            
            # calculate change in energy associated with spin flip
            ΔEᵢⱼ = calc_ΔEᵢⱼ(s,i,j,J,B)
            
            # calculate metropolis-hastings acceptance probability
            Pᵢⱼ  = min(1.0, exp(-ΔEᵢⱼ/T))
            
            # accept or reject proposed spin flip
            if rand() < Pᵢⱼ
                s[i,j] = -s[i,j]
            end
        end
    end
    
    return nothing
end

function run_simulation()

    # system size
    L = 100

    # temperature
    T = 2.5

    # coupling strength
    J = 1.0

    # initial magnetic field
    B = 0.0

    # number of Monte Carlo sweeps to perform
    N = 100_000

    # target magnetization
    m₀ = 0.5

    # initialize ising spins to random initial configuration
    s = rand(-1:2:1,L,L)

    # initialize MuTunerLogger for tuning magnetic field B
    Btuner = MuTunerLogger(n₀=m₀, β=1/T, V=L^2, u₀=J, μ₀=B, c=0.5, s=1.0+0.0im)

    # perform Monte Carlo simulation
    for i in 1:N
        
        # perform Monte Carlo sweep through lattice
        monte_carlo_sweep!(s, T, J, B)
        
        # measure ⟨m⟩ = ⟨M⟩/L²
        M = float(sum(s))
        m = M/L^2
        
        # measure ⟨M²⟩
        M² = M^2
        
        # update the magnetic field
        B = update!(μtuner=Btuner, n=m, N²=M², s=1.0+0.0im)
    end

    # replay the tuning time-series for all relevant observables and write to file.
    save(Btuner, "Btuner.csv")

    return nothing
end

# run the simulation
run_simulation()
```