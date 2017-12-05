#==========================================================================================#

# INFLUENCE FUNCTION

function influence(obj::ParModel, w::AbstractWeights)
    s = score(obj)
    j = jacobian(obj, w)
    return scale!(- 1.0, s / j')
end

function influence(obj::TwoStageModel, w::AbstractWeights)
    c = influence(first_stage(obj), w) * crossjacobian(obj, w)'
    s = score(obj)
    j = jacobian(obj, w)
    return scale!(- 1.0, (s + c) / j')
end

#==========================================================================================#

# COVARIANCE MATRIX FOR HOMOSCEDASTIC MLE

function _vcov!(obj::MLE, corr::Homoscedastic, w::AbstractWeights)
    if getcorr(obj).method == "OIM"
        obj.V = scale!(- 1.0, inv(jacobian(obj, w)))
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
    Ψ     = scale!(w, copy(ψ))
    obj.V = ψ' * Ψ
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::Heteroscedastic, w::AbstractWeights)
    ψ     = influence(obj, w) ; scale!(w, ψ)
    obj.V = crossprod(ψ)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::ClusterOrCross, w::UnitWeights)
    ψ     = influence(obj, w)
    obj.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.V, obj, corr)
end

function _vcov!(obj::ParModel, corr::ClusterOrCross, w::AbstractWeights)
    ψ     = influence(obj, w) ; scale!(w, ψ)
    obj.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.V, obj, corr)
end

# COVARIANCE MATRIX FOR TWO-STAGE PARAMETRIC MODELS

function _vcov!(obj::TwoStageModel, corr::Heteroscedastic, w::UnitWeights)
    ψ                  = influence(obj, w)
    obj.second_stage.V = crossprod(ψ)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::Heteroscedastic, w::FrequencyWeights)
    ψ                  = influence(obj, w)
    Ψ                  = scale!(w, copy(ψ))
    obj.second_stage.V = ψ' * Ψ
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::Heteroscedastic, w::AbstractWeights)
    ψ                  = influence(obj, w) ; scale!(w, ψ)
    obj.second_stage.V = crossprod(ψ)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::ClusterOrCross, w::UnitWeights)
    ψ                  = influence(obj, w)
    obj.second_stage.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.second_stage.V, obj, corr)
end

function _vcov!(obj::TwoStageModel, corr::ClusterOrCross, w::AbstractWeights)
    ψ                  = influence(obj, w) ; scale!(w, ψ)
    obj.second_stage.V = crossprod(ψ, corr.mat)
    adjfactor!(obj.second_stage.V, obj, corr)
end
