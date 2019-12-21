#==========================================================================================#

# INFLUENCE FUNCTION

function influence(obj::MLE, w::AbstractWeights)
    s = score(obj)
    j = jacobian(obj, w)
    return - s / j'
end

function influence(obj::TwoStageModel, w::AbstractWeights)
    c = influence(first_stage(obj), w) * crossjacobian(obj, w)'
    s = score(obj)
    j = jacobian(obj, w)
    return - (s + c) / j'
end

#==========================================================================================#

# FINITE-SAMPLE ADJUSTMENT

varcorrection(corr::Homoscedastic, w::AbstractWeights) = 1

function varcorrection(corr::Heteroscedastic, w::AbstractWeights)
    sum(w) * varcorrection(w, true)
end

varcorrection(corr::Clustered, w::AbstractWeights) = corr.n_clusters / (corr.n_clusters - 1)

function varcorrection(corr::CrossCorrelated, w::AbstractWeights)
    sum(w) * varcorrection(w, true)
end

# SET COVARIANCE MATRIX

_vcov!(obj::ParModel, V::AbstractMatrix) = (switch_stage(obj).V = V)

#==========================================================================================#

# COVARIANCE MATRIX FOR HOMOSCEDASTIC MLE

function _vcov!(obj::MLE, corr::Homoscedastic, w::AbstractWeights)
    obj.V = getcorr(obj).expected ? inv(crossprod(score(obj), w)) : - inv(jacobian(obj, w))
end

#==========================================================================================#

# COVARIANCE MATRIX FOR MLE AND TWO-STAGE MODELS

for m in (MLE, TwoStageModel)
    @eval begin

        function _vcov!(obj::$m, corr::Heteroscedastic, w::UnitWeights)
            ψ = influence(obj, w)
            V = crossprod(ψ)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Heteroscedastic, w::FrequencyWeights)
            ψ = influence(obj, w)
            V = crossprod(ψ, w)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Heteroscedastic, w::AbstractWeights)
            ψ = influence(obj, w) ; lmul!(Diagonal(w), ψ)
            V = crossprod(ψ)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Clustered, w::UnitWeights)
            ψ = corr.matrix * influence(obj, w)
            V = crossprod(ψ)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Clustered, w::AbstractWeights)
            ψ = corr.matrix * influence(obj, w) ; lmul!(Diagonal(w), ψ)
            V = crossprod(ψ)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::CrossCorrelated, w::UnitWeights)
            ψ = influence(obj, w)
            V = crossprod(ψ, corr.matrix)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::CrossCorrelated, w::AbstractWeights)
            ψ = influence(obj, w) ; lmul!(Diagonal(w), ψ)
            V = crossprod(ψ, corr.matrix)
            corr.corrected && lmul!(varcorrection(corr, w), V)
            _vcov!(obj, V)
        end
    end
end

#==========================================================================================#

# GMM WEIGHT MATRIX

wmatrix(obj::GMM, corr::Heteroscedastic, ::UnitWeights) = crossprod(score(obj))
wmatrix(obj::GMM, corr::Clustered, ::UnitWeights) = crossprod(corr.matrix * score(obj))
wmatrix(obj::GMM, corr::CrossCorrelated, ::UnitWeights) = crossprod(score(obj), corr.matrix)
wmatrix(obj::GMM, corr::Heteroscedastic, w::FrequencyWeights) = crossprod(score(obj), w)

function wmatrix(obj::GMM, corr::Heteroscedastic, w::AbstractWeights)
    s = score(obj) ; lmul!(Diagonal(w), s)
    crossprod(s)
end

function wmatrix(obj::GMM, corr::Clustered, w::AbstractWeights)
    s = corr.matrix * score(obj) ; lmul!(Diagonal(w), s)
    crossprod(s)
end

function wmatrix(obj::GMM, corr::CrossCorrelated, w::AbstractWeights)
    s = score(obj) ; lmul!(Diagonal(w), s)
    crossprod(s, corr.matrix)
end

# COVARIANCE MATRIX FOR GMM

function _vcov!(obj::GMM, corr::CorrStructure, w::AbstractWeights)

    j = jacobian(obj, w)
    S = wmatrix(obj, corr, w)

    if obj.method == "Method of moments"
        V = j \ (S / j')
    elseif obj.method == "One-step GMM"
        ψ = j / (j' * j)
        V = crossprod(ψ, S)
    elseif obj.method == "Two-step GMM"
        ω = obj.W \ j
        ψ = ω / (j' * ω)
        V = crossprod(ψ, S)
    elseif obj.method == "Optimal GMM"
        ω = S \ j
        V = inv(j' * ω)
    end

    corr.corrected && lmul!(varcorrection(corr, w), V)
    _vcov!(obj, V)
end
