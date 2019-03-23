#==========================================================================================#

# This file replicates tests from the official documentation of Stata
# Some DOF adjustements are necessary because Stata uses different corrections

#==========================================================================================#

# MISE EN PLACE

using Test
using CSV
using DataFrames
using Microeconometrics
using StatsBase
using StatsModels

const datadir = joinpath(dirname(@__FILE__), "..", "data")

#==========================================================================================#

S        = CSV.read(joinpath(datadir, "auto.csv"))
S[:gpmw] = ((1.0 ./ S[:mpg]) ./ S[:weight]) * 100000
M        = Dict(:response => "gpmw", :control => "foreign + 1")

@testset "OLS Homoscedastic" begin

    D = Microdata(S, M, vcov = Homoscedastic())
    E = fit(OLS, D)

    β = [0.24615258; 1.60900394]
    σ = [0.05494872; 0.02996078]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(loglikelihood(E), 9.39839211)
    @test isapprox(nullloglikelihood(E), 0.30171734)
    @test isapprox(deviance(E), 3.36079458)
    @test isapprox(nulldeviance(E), 4.29750019)
    @test isapprox(r2(E), 0.21796522, atol = 1e-7, rtol = 1e-7)
    @test isapprox(adjr2(E), 0.20710363, atol = 1e-7, rtol = 1e-7)
    @test dof(E) == 2
end

@testset "OLS Heterosc." begin

    D = Microdata(S, M, vcov = Heteroscedastic())
    E = fit(OLS, D)

    β = [0.24615258; 1.60900394]
    σ = [0.06792384; 0.02345347]
    c = (nobs(E) - dof(E)) / (nobs(E) - 1)
    σ = σ * sqrt(c)

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(r2(E), 0.21796522, atol = 1e-7, rtol = 1e-7)
    @test isapprox(adjr2(E), 0.20710363, atol = 1e-7, rtol = 1e-7)
end

T        = Dict("age" => Union{Int, Missing})
S        = CSV.read(joinpath(datadir, "regsmpl.csv"), types = T)
S[:age2] = Array{eltype(S[:age])}(S[:age].^2)
W        = Clustered(S, :idcode)
M        = Dict(:response => "ln_wage", :control => "age + age2 + tenure + 1")

@testset "OLS Clustered" begin

    D = Microdata(S, M, vcov = W)
    E = fit(OLS, D)

    β = [0.07521723; -0.00108513; 0.03908767; 0.33398213]
    σ = [0.00457113;  0.00007785; 0.00144254; 0.06419184]
    c = (nobs(E) - dof(E)) / (nobs(E) - 1)
    σ = σ * sqrt(c)

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(r2(E), 0.16438516, atol = 1e-7, rtol = 1e-7)
    @test isapprox(adjr2(E), 0.16429594, atol = 1e-7, rtol = 1e-7)
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

    D = Microdata(S, M, vcov = Homoscedastic())
    E = fit(IV, D)

    β = [0.00223983; 0.08151597; 120.70651454]
    σ = [0.00033876; 0.30815277;  15.70688390]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(r2(E), 0.59888202, atol = 1e-7, rtol = 1e-7)
    @test isapprox(adjr2(E), 0.58181317, atol = 1e-7, rtol = 1e-7)
end

@testset "IV GMM" begin

    D = Microdata(S, M, vcov = Heteroscedastic()) ;
    E = fit(IV, D, method = "Two-step GMM")

    β = [0.00146433; 0.76154816; 112.12271295] ;
    σ = [0.00044727; 0.28951046;  10.80234023] ;
    c = nobs(E) / (nobs(E) - 1) ;
    σ = σ * sqrt(c) ;

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(r2(E), 0.66162834, atol = 1e-7, rtol = 1e-7)
    @test isapprox(adjr2(E), 0.64722954, atol = 1e-7, rtol = 1e-7)
end

#==========================================================================================#

S = CSV.read(joinpath(datadir, "lbw.csv"))
M = Dict(:response => "low", :control => "age + lwt + race + smoke + ptl + ht + ui + 1")
C = Dict(:race => DummyCoding(base = "white"))
D = Microdata(S, M, vcov = Homoscedastic(), contrasts = C)

@testset "Logit" begin

    E = fit(Logit, D)

    β = [-0.02710031; -0.01515082; 1.26264728; 0.86207916; 0.92334482;
          0.54183656;  1.83251780; 0.75851348; 0.46122388]
    σ = [ 0.03645043;  0.00692588; 0.52641014; 0.43915315; 0.40082664;
          0.34624900; 0.69162923; 0.45937677; 1.20458975]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(deviance(E), 201.44799113, atol = 1e-7, rtol = 1e-7)
    @test isapprox(loglikelihood(E), -100.72399557, atol = 1e-7, rtol = 1e-7)
    @test isapprox(nullloglikelihood(E), -117.33599810, atol = 1e-7, rtol = 1e-7)
    @test isapprox(aic(E), 219.44799113, atol = 1e-7, rtol = 1e-7)
    @test isapprox(aicc(E), 220.45357772, atol = 1e-7, rtol = 1e-7)
    @test isapprox(bic(E), 248.62371427, atol = 1e-7, rtol = 1e-7)
end

@testset "Probit" begin

    E = fit(Probit, D)

    β = [-0.01754447; -0.00882045; 0.74752563; 0.51447107; 0.56276006;
          0.31782665;  1.09945075; 0.46279438; 0.26827531]
    σ = [ 0.02162924;  0.00397343; 0.31664054; 0.25558235; 0.23577825;
          0.20012526;  0.41927933; 0.27560931; 0.70152540]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(deviance(E), 201.12189887, atol = 1e-7, rtol = 1e-7)
    @test isapprox(loglikelihood(E), -100.56094943, atol = 1e-7, rtol = 1e-7)
    @test isapprox(nullloglikelihood(E), -117.33599810, atol = 1e-7, rtol = 1e-7)
end

@testset "Cloglog" begin

    E = fit(Cloglog, D)

    β = [-0.02309612; -0.01129839; 1.07937099; 0.72856146; 0.73324598;
          0.33131635;  1.42557598; 0.56457925; -0.09224780]
    σ = [0.02906676; 0.00521572; 0.40185917; 0.32942310; 0.30247340;
         0.20833347; 0.45334563; 0.34357409; 0.90994949]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(deviance(E), 202.16661013, atol = 1e-7, rtol = 1e-7)
    @test isapprox(loglikelihood(E), -101.08330506, atol = 1e-7, rtol = 1e-7)
    @test isapprox(nullloglikelihood(E), -117.33599810, atol = 1e-7, rtol = 1e-7)
end

@testset "Gompit" begin

    E = fit(Gompit, D)

    β = [ 0.01937153; 0.00827072; -0.65246666; -0.46606001; -0.55757061;
         -0.40377232; -1.06494339; -0.51695211; -0.63626157]
    σ = [0.02104003; 0.00375728; 0.32497917; 0.24843060; 0.23366182;
         0.25411648; 0.48787151; 0.29917807; 0.67163691]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(deviance(E), 200.35974543, atol = 1e-7, rtol = 1e-7)
    @test isapprox(loglikelihood(E), -100.17987271, atol = 1e-7, rtol = 1e-7)
    @test isapprox(nullloglikelihood(E), -117.33599810, atol = 1e-7, rtol = 1e-7)
end

#==========================================================================================#

S          = CSV.read(joinpath(datadir, "dollhill3.csv"))
S[:pyears] = log.(S[:pyears])

M = Dict(
        :response => "deaths",
        :control  => "smokes + agecat + 1",
        :offset   => "pyears"
    )

@testset "Poisson" begin

    D = Microdata(S, M, vcov = Homoscedastic())
    E = fit(Poisson, D)

    β = [0.35453564; 1.48400701; 2.62750512; 3.35049279; 3.70009645; -7.91932571]
    σ = [0.10737412; 0.19510337; 0.18372727; 0.18479918; 0.19221951;  0.19176182]

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test isapprox(deviance(E), 12.13236640, atol = 1e-7, rtol = 1e-7)
    @test isapprox(nulldeviance(E), 644.26899081, atol = 1e-7, rtol = 1e-7)
    @test isapprox(loglikelihood(E), -33.60015344, atol = 1e-7, rtol = 1e-7)
    @test isapprox(nullloglikelihood(E), -495.06763568, atol = 1e-7, rtol = 1e-7)
end

S = CSV.read(joinpath(datadir, "website.csv"))
M = Dict(
        :response   => "visits",
        :control    => "ad + female + 1",
        :treatment  => "time",
        :instrument => "phone + frfam"
    )

@testset "IVPoisson" begin

    D = Microdata(S, M)
    E = fit(IVPoisson, D, method = "Two-step GMM")

    β = [0.05892941; 0.13734403; -0.02477071; 1.04150546]
    σ = [0.01079421; 0.01015705;  0.03762177; 0.03858479]
    c = nobs(E) / (nobs(E) - 1)
    σ = σ * sqrt(c)

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
end

S = CSV.read(joinpath(datadir, "trip.csv"))
M = Dict(
        :response   => "trips",
        :control    => "cbd + ptn + worker + weekend + 1",
        :treatment  => "tcost",
        :instrument => "pt"
    )

@testset "Mullahy" begin

    D = Microdata(S, M)
    E = fit(Mullahy, D, method = "Two-step GMM")

    β = [0.03521846; -0.00839804; -0.01131463; 0.66230176; 0.30093231; 0.26544230]
    σ = [0.00981818;  0.00201721;  0.00218185; 0.05199088; 0.03626820; 0.15501267]
    c = nobs(E) / (nobs(E) - 1)
    σ = σ * sqrt(c)

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
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
    σ = [  25.81524380;    9.57136890]
    c = nobs(E) / (nobs(E) - 1)
    σ = σ * sqrt(c)

    @test isapprox(coef(E), β, atol = 1e-7, rtol = 1e-7)
    @test isapprox(stderror(E), σ, atol = 1e-7, rtol = 1e-7)
    @test dof(E) == 2
end

#==========================================================================================#

S = CSV.read(joinpath(datadir, "income.csv"))
X = (S[:male] .== 1)
M = Dict(:response => "inc", :control => "edu + exp + 1")

@testset "Hausman" begin

    E0    = fit(OLS, Microdata(S, M, subset = X))
    E0.V .= E0.V .* ((nobs(E0) - 1) / nobs(E0)) * (277 / 276)
    E1    = fit(OLS, Microdata(S, M, subset = .!X))
    E1.V .= E1.V .* ((nobs(E1) - 1) / nobs(E1)) * (277 / 276)
    E     = hausman_2s(E0, E1)

    p = [0.20345785; 0.59592600; 0.28068223]

    @test isapprox(pval(E), p, atol = 1e-7, rtol = 1e-7)
end
