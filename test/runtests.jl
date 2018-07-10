#==========================================================================================#

# This file replicates tests from the official documentation of Stata
# Some DOF adjustements are necessary because Stata uses different corrections

#==========================================================================================#

# MISE EN PLACE

using Base.Test
using CSV
using DataFrames
using Microeconometrics

function test_show(x)
    io = IOBuffer()
    show(io, x)
end

const datadir = joinpath(dirname(@__FILE__), "..", "data")

#==========================================================================================#

S        = CSV.read(joinpath(datadir, "auto.csv"))
S[:gpmw] = ((1.0 ./ S[:mpg]) ./ S[:weight]) * 100000
M        = Dict(:response => "gpmw", :control => "foreign + 1")

@testset "OLS_Homoscedastic" begin

    D = Microdata(S, M, vcov = Homoscedastic())
    E = fit(OLS, D)

    β = [0.24615258;  1.60900394]
    t = [4.47967818; 53.70367955]

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test isapprox(r2(E), 0.21796522, atol = 1e-7)
    @test isapprox(adjr2(E), 0.20710363, atol = 1e-7)
    @test dof(E) == 2
end

@testset "OLS_Heteroscedastic" begin
    
    D = Microdata(S, M, vcov = Heteroscedastic())
    E = fit(OLS, D)

    β = [0.24615258;  1.60900394]
    t = [3.62395001; 68.60409513]
    c = (nobs(D) - dof(E)) / (nobs(D) - 1)
    t = t / sqrt(c)

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test isapprox(r2(E), 0.21796522, atol = 1e-7)
    @test isapprox(adjr2(E), 0.20710363, atol = 1e-7)
    @test dof(E) == 2
end

T        = Dict("age" => Union{Int, Missing})
S        = CSV.read(joinpath(datadir, "regsmpl.csv"), types = T)
S[:age2] = Array{eltype(S[:age])}(S[:age].^2)
W        = Clustered(S, :idcode)
M        = Dict(:response => "ln_wage", :control => "age + age2 + tenure + 1")

@testset "OLS_Clustered" begin

    D = Microdata(S, M, vcov = W)
    E = fit(OLS, D)

    β = [ 0.07521723; - 0.00108513;  0.03908767; 0.33398213]
    t = [16.45483785; -13.93935551; 27.09636359; 5.20287523]
    c = (nobs(D) - dof(E)) / (nobs(D) - 1)
    t = t / sqrt(c)

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test isapprox(r2(E), 0.16438516, atol = 1e-7)
    @test isapprox(adjr2(E), 0.16429594, atol = 1e-7)
    @test dof(E) == 4
end

#==========================================================================================#

S = CSV.read(joinpath(datadir, "lbw.csv"))
M = Dict(:response => "low", :control => "age + lwt + race + smoke + ptl + ht + ui + 1")
C = Dict(:race => DummyCoding(base = "white"))

@testset "Logit" begin

    D = Microdata(S, M, vcov = Homoscedastic(), contrasts = C)
    E = fit(Logit, D)

    β = [-0.02710031; -0.01515082; 1.26264728; 0.86207916; 0.92334482;
          0.54183656;  1.83251780; 0.75851348; 0.46122388]
    t = [-0.74348404; -2.18756627; 2.39859983; 1.96304899; 2.30360141;
          1.56487547;  2.64956671; 1.65117944; 0.38288876]

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test isapprox(deviance(E), 201.44799113, atol = 1e-7)
    @test isapprox(loglikelihood(E), -100.72399557, atol = 1e-7)
    @test isapprox(aic(E), 219.44799113, atol = 1e-7)
    @test isapprox(aicc(E), 220.45357772, atol = 1e-7)
    @test isapprox(bic(E), 248.62371427, atol = 1e-7)
    @test dof(E) == 9
end

@testset "Probit" begin

    E = fit(Probit, D) ;

    β = [-0.01754447; -0.00882045; 0.74752563; 0.51447107; 0.56276006;
          0.31782665; 1.09945075; 0.46279438; 0.26827531]
    t = [-0.81114596; -2.21985717; 2.36080203; 2.01293659; 2.38681919;
          1.58813858; 2.62223933; 1.67916817; 0.38241710]

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test isapprox(deviance(E), 201.12189887, atol = 1e-7)
    @test isapprox(loglikelihood(E), -100.56094943, atol = 1e-7)
    @test isapprox(aic(E), 219.12189887, atol = 1e-7)
    @test isapprox(aicc(E), 220.12748546, atol = 1e-7)
    @test isapprox(bic(E), 248.29762201, atol = 1e-7)
    @test dof(E) == 9
end

#==========================================================================================#

S         = CSV.read(joinpath(datadir, "cattaneo2.csv"))
S[:mage2] = S[:mage].^2

M = Dict(
        :response  => "bweight",
        :control   => "mmarried + mage + mage2 + fbaby + medu + 1",
        :treatment => "mbsmoke",
    )

@testset "IPW" begin

    D = Microdata(S, M)
    E = fit(IPW, Probit, D)

    β = [-230.68863780; 3403.46270868]
    t = [  -8.93614020; 355.58787307]
    c = nobs(E) / (nobs(E) - 1)
    t = t / sqrt(c)

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test dof(E) == 2
end

#==========================================================================================#

S = CSV.read(joinpath(datadir, "hsng2.csv"))
M = Dict(
        :response   => "rent",
        :control    => "pcturban + 1",
        :treatment  => "hsngval",
        :instrument => "faminc + region"
    )

@testset "TSLS" begin

    D = Microdata(S, M)
    E = fit(IV, D)

    β = [0.00223983; 0.08151597]
    t = [3.33306931; 0.18334930]
    c = nobs(E) / (nobs(E) - 1)
    t = t / sqrt(c)

    test_show(E)
    @test isapprox(coef(E), β, atol = 1e-7)
    @test isapprox(tstat(E), t, atol = 1e-7)
    @test isapprox(r2(E), 0.59888202, atol = 1e-7)
    @test isapprox(adjr2(E), 0.58181317, atol = 1e-7)
    @test dof(E) == 2
end
