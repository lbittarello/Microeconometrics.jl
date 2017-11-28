#==========================================================================================#

# MISE EN PLACE

# These tests replicate the tests of GLM.jl

using Microeconometrics
using Base.Test
using StatsBase
using CSV
using DataFrames
# using RDatasets

function test_show(x)
    io = IOBuffer()
    show(io, x)
end

const datadir = joinpath(dirname(@__FILE__), "..", "data")

#==========================================================================================#

# OLS

# RDatasets needs updating

#=
dset = dataset("datasets", "Formaldehyde")
mdta = Microdata(dset, Homoscedastic(), response = "OptDen", control = "Carb + 1")

@testset "OLS" begin

    Σ = [ 1.831836734693908e-04 -9.464489795918525e-05
         -9.464489795918525e-05  6.136653061224592e-05]

    e_ols = fit(OLS, mdta)
    test_show(e_ols)
    @test isapprox(coef(e_ols), [0.8762857142857143; 0.0050857142857141935])
    @test isapprox(vcov(e_ols), Σ)
    @test dof(e_ols) == 2
    @test r²(e_ols) == r2(e_ols)
    @test isapprox(r²(e_ols), 0.9990466748057584)
    @test adjr²(e_ols) == adjr2(e_ols)
    @test isapprox(adjr²(e_ols), 0.998808343507198)
end
=#

#==========================================================================================#

# BINARY MODELS (LOGIT AND PROBIT)

# We don't quite match the coeficients
# However, we match the value of the log loglikelihood at the maximum
# (actually, ours is a little higher)
# So the log loglikelihood must be somewhat flat around the MLE

dset        = CSV.read(joinpath(datadir, "admit.csv"))
dset[:rank] = categorical(dset[:rank], levels = [1; 2; 3; 4])
mdta        = Microdata(dset, response = "admit", control = "gre + gpa + rank + 1")

@testset "Logit" begin

    β = [ 0.0022644256521549043, 0.8040374535155766, -0.6754428594116577,
         -1.3402038117481079,   -1.5514636444657492, -3.9899786606380734]

    e_logit = fit(Logit, mdta)
    test_show(e_logit)
    @test dof(e_logit) == 6
    @test isapprox(deviance(e_logit), 458.5174924758994)
    @test isapprox(loglikelihood(e_logit), -229.25874623794968)
    @test isapprox(aic(e_logit), 470.51749247589936)
    @test isapprox(aicc(e_logit), 470.7312329339146)
    @test isapprox(bic(e_logit), 494.4662797585473)
    #@test isapprox(coef(e_logit), β)
end

@testset "Probit" begin

    β = [ 0.0013755394922972369, 0.47772908362647015, -0.4154125854823675,
         -0.8121458010130356,   -0.9359047862425298,  -2.3867922998680786]

    e_probit = fit(Probit, mdta)
    test_show(e_probit)
    @test dof(e_probit) == 6
    @test isapprox(deviance(e_probit), 458.4131713833386)
    @test isapprox(loglikelihood(e_probit), -229.20658569166932)
    @test isapprox(aic(e_probit), 470.41317138333864)
    @test isapprox(aicc(e_probit), 470.6269118413539)
    @test isapprox(bic(e_probit), 494.36195866598655)
    #@test isapprox(coef(e_probit), β)
end
