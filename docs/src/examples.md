# Examples

On this page we include examples for how this package can be used.

## The Ising Model

While the algorithm implmeneted in this package, introduced in the paper
[arXiv:2201.01296](https://arxiv.org/abs/2201.01296), was initially designed to tune the
chemical potential ``\mu`` to achieve a target particle number ``\langle N \rangle``, it can
also be applied to achieve a target magnetization in a spin model by tuning the
applied magnetic field.

In this example we use [`MuTuner.jl`](https://github.com/cohensbw/MuTuner.jl) to tune the
magnetic field ``B`` to achieve a target magnetization ``\langle m \rangle`` in a Monte Carlo
simulation of the Ising model
```math
H = -J \sum_{\langle ij\rangle} s_{i}s_{j} - B \sum_{i} s_{i},
```
where ``J`` is the coupling strength.
The tuning algorithm can be directly applied by applying the substitutions
```math
\mu \mapsto B,\quad N \mapsto M,\quad \kappa \mapsto \chi,
```
where ``\kappa`` and ``\chi`` are the compressibility and magnetic
susceptibility respectively.

Below is an example script using [`MuTuner.jl`](https://github.com/cohensbw/MuTuner.jl) in a
Monte Carlo simulation of the Ising Model:

```julia
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
    for j in 1:L
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
Btuner = MuTunerLogger(n₀=m₀, β=1/T, V=L^2, U₀=J, μ₀=B, c=0.5)

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
    B = update!(μtuner=Btuner, n=m, N²=M²)
end

# replay the tuning timeseries for all relevant observables.
# these resulting timeseries arrays can be used to generate a figure
# comprable to Figure 1 in the paper arXiv:2201.01296.
tape        = replay(Btuner);
B_traj      = tape.μ_traj
B_bar_traj  = tape.μ_bar_traj
m_traj      = tape.N_traj / L^2
m_bar_traj  = tape.N_bar_traj / L^2
χ_bar_traj  = tape.κ_bar_traj
χ_fluc_traj = tape.κ_fluc_traj
χ_min_traj  = tape.κ_min_traj
χ_max_traj  = tape.κ_max_traj
```