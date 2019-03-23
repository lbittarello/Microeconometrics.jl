#==========================================================================================#

# INFLUENCE FUNCTION

function influence(obj::ParModel, w::AbstractWeights)
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

function influence(obj::GMM, w::AbstractWeights)

    s = score(obj)
    j = jacobian(obj, w)

    if obj.method == "Method of moments"
        return - s / j'
    elseif obj.method == "One-step GMM"
        ψ = j / (j' * j)
        return - s * ψ
    elseif (obj.method == "Two-step GMM") | (obj.method == "Optimal GMM")
        ω = obj.W \ j
        ψ = ω / (j' * ω)
        return - s * ψ
    end
end

#==========================================================================================#

# SET COVARIANCE MATRIX

_vcov!(obj::ParModel, V::Matrix)      = (obj.V = V)
_vcov!(obj::GMM, V::Matrix)           = (obj.V = V)
_vcov!(obj::TwoStageModel, V::Matrix) = (obj.second_stage.V = V)

#==========================================================================================#

# COVARIANCE MATRIX FOR HOMOSCEDASTIC MLE

function _vcov!(obj::MLE, corr::Homoscedastic, w::AbstractWeights)
    if getcorr(obj).method == "OIM"
        obj.V = - inv(jacobian(obj, w))
    elseif getcorr(obj).method == "OPG"
        g     = score(obj)
        obj.V = inv(crossprod(g, w))
    else
        throw("unknown method")
    end
end

#==========================================================================================#

# COVARIANCE MATRIX FOR MLE AND TWO-STAGE MODELS

for m in (MLE, TwoStageModel)
    @eval begin

        function _vcov!(obj::$m, corr::Heteroscedastic, w::UnitWeights)
            ψ = influence(obj, w)
            V = crossprod(ψ)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Heteroscedastic, w::FrequencyWeights)
            ψ = influence(obj, w)
            V = crossprod(ψ, w)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Heteroscedastic, w::AbstractWeights)
            ψ = influence(obj, w) ; lmul!(Diagonal(w), ψ)
            V = crossprod(ψ)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Clustered, w::UnitWeights)
            ψ = corr.mat * influence(obj, w)
            V = crossprod(ψ)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::Clustered, w::AbstractWeights)
            ψ = corr.mat * influence(obj, w) ; lmul!(Diagonal(w), ψ)
            V = crossprod(ψ)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::CrossCorrelated, w::UnitWeights)
            ψ = influence(obj, w)
            V = crossprod(ψ, corr.mat)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end

        function _vcov!(obj::$m, corr::CrossCorrelated, w::AbstractWeights)
            ψ = influence(obj, w) ; lmul!(Diagonal(w), ψ)
            V = crossprod(ψ, corr.mat)
            adjfactor!(V, obj, corr)
            _vcov!(obj, V)
        end
    end
end

#==========================================================================================#

# GMM WEIGHT MATRIX

function wmatrix(obj::GMM, corr::Heteroscedastic, ::UnitWeights)
    s = score(obj)
    return crossprod(s)
end

function wmatrix(obj::GMM, corr::Heteroscedastic, w::FrequencyWeights)
    s = score(obj)
    return crossprod(s, w)
end

function wmatrix(obj::GMM, corr::Heteroscedastic, w::AbstractWeights)
    s = score(obj) ; lmul!(Diagonal(w), s)
    return crossprod(s)
end

function wmatrix(obj::GMM, corr::Clustered, ::UnitWeights)
    s = corr.mat * score(obj)
    return crossprod(s)
end

function wmatrix(obj::GMM, corr::Clustered, w::AbstractWeights)
    s = corr.mat * score(obj) ; lmul!(Diagonal(w), s)
    return crossprod(s)
end

function wmatrix(obj::GMM, corr::CrossCorrelated, ::UnitWeights)
    s = score(obj)
    return crossprod(s, corr.mat)
end

function wmatrix(obj::GMM, corr::CrossCorrelated, w::AbstractWeights)
    s = score(obj) ; lmul!(Diagonal(w), s)
    return crossprod(s, corr.mat)
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

    adjfactor!(V, obj, corr)
    _vcov!(obj, V)
end
