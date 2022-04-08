using MuTuner
using Test
using Statistics

@testset "MuTuner.jl" begin
    
    # test forgetful mean and variance
    c      = 0.5
    N      = 10
    x      = rand(N)
    x̄_Ot   = 0.0
    v̄_Ot   = 0.0
    x̄_Ot²  = 0.0
    v̄_Ot²  = 0.0
    x̄s_Ot  = zeros(N)
    v̄s_Ot  = zeros(N)
    x̄s_Ot² = zeros(N)
    v̄s_Ot² = zeros(N)
    for t in 0:(N-1)
        # update forgetful mean and variance
        xs_Ot       = @view x[1:t+1]
        x̄_Ot, v̄_Ot  = update_forgetful_mv(xs_Ot, x̄_Ot, v̄_Ot, c)
        x̄s_Ot[t+1]  = x̄_Ot
        v̄s_Ot[t+1]  = v̄_Ot
        # calculate forgetful mean and variance from scratch
        xs_Ot²      = @view x[ceil(Int,c*t)+1:t+1]
        x̄_Ot²       = mean(xs_Ot²)
        v̄_Ot²       = var(xs_Ot², corrected=false)
        x̄s_Ot²[t+1] = x̄_Ot²
        v̄s_Ot²[t+1] = v̄_Ot²
    end
    # test that forgetful mean is calculated correctly
    @test all(x̄s_Ot .≈ x̄s_Ot²)
    # test that forgetful variance is calculated correctly
    @test all(v̄s_Ot .≈ v̄s_Ot²)

end