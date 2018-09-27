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

# COVARIANCE MATRIX FOR ONE-STAGE PARAMETRIC MODELS

function _vcov!(obj::ParModel, corr::Heteroscedastic, w::UnitWeights)
    ψ     = influence(obj, w)
    obj.V = crossprod(ψ)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::Heteroscedastic, w::FrequencyWeights)
    ψ     = influence(obj, w)
    obj.V = crossprod(ψ, w)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::Heteroscedastic, w::AbstractWeights)
    ψ     = influence(obj, w) ; lmul!(Diagonal(w), ψ)
    obj.V = crossprod(ψ)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::Clustered, w::UnitWeights)
    ψ     = corr.mat * influence(obj, w)
    obj.V = crossprod(ψ)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::Clustered, w::AbstractWeights)
    ψ     = corr.mat * influence(obj, w) ; lmul!(Diagonal(w), ψ)
    obj.V = crossprod(ψ)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::CrossCorrelated, w::UnitWeights)
    ψ     = influence(obj, w)
    obj.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::CrossCorrelated, w::AbstractWeights)
    ψ     = influence(obj, w) ; lmul!(Diagonal(w), ψ)
    obj.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.V, obj, corr)
end

#==========================================================================================#

# COVARIANCE MATRIX FOR TWO-STAGE PARAMETRIC MODELS

function _vcov!(obj::TwoStageModel, corr::Heteroscedastic, w::UnitWeights)
    ψ                  = influence(obj, w)
    obj.second_stage.V = crossprod(ψ)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::Heteroscedastic, w::FrequencyWeights)
    ψ                  = influence(obj, w)
    obj.second_stage.V = crossprod(ψ, w)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::Heteroscedastic, w::AbstractWeights)
    ψ                  = influence(obj, w) ; lmul!(Diagonal(w), ψ)
    obj.second_stage.V = crossprod(ψ)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::Clustered, w::UnitWeights)
    ψ                  = corr.mat * influence(obj, w)
    obj.second_stage.V = crossprod(ψ)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::Clustered, w::AbstractWeights)
    ψ                  = corr.mat * influence(obj, w) ; lmul!(Diagonal(w), ψ)
    obj.second_stage.V = crossprod(ψ)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::CrossCorrelated, w::UnitWeights)
    ψ                  = influence(obj, w)
    obj.second_stage.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::CrossCorrelated, w::AbstractWeights)
    ψ                  = influence(obj, w) ; lmul!(Diagonal(w), ψ)
    obj.second_stage.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.second_stage.V, obj, corr)
end

#==========================================================================================#

# GMM WEIGHT MATRIX

function wmatrix(obj::GMM, corr::Heteroscedastic, w::UnitWeights)
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

function wmatrix(obj::GMM, corr::Clustered, w::UnitWeights)
    s = corr.mat * score(obj)
    return crossprod(s)
end

function wmatrix(obj::GMM, corr::Clustered, w::AbstractWeights)
    s = corr.mat * score(obj) ; lmul!(Diagonal(w), s)
    return crossprod(s)
end

function wmatrix(obj::GMM, corr::CrossCorrelated, w::UnitWeights)
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
        obj.V = j \ (S / j')
    elseif obj.method == "One-step GMM"
        ψ     = j / (j' * j)
        obj.V = crossprod(ψ, S)
    elseif obj.method == "Two-step GMM"
        ω     = obj.W \ j
        ψ     = ω / (j' * ω)
        obj.V = crossprod(ψ, S)
    elseif obj.method == "Optimal GMM"
        ω     = S \ j
        obj.V = inv(j' * ω)
    end

    adjfactor!(obj.V, obj, corr)
end
