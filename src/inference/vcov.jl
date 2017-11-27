
#==========================================================================================#

# INFLUENCE FUNCTION

function influence(obj::ParModel, args...)
    s = score(obj, args...)
    j = jacobian(obj, args...)
    return scale!(- 1.0, s / j')
end

function influence(obj::TwoStageModel, args...)
    c = influence(first_stage(obj), args...) * crossjacobian(obj, args...)'
    s = score(obj, args...)
    j = jacobian(obj, args...)
    return scale!(- 1.0, (s + c) / j')
end

#==========================================================================================#

# COVARIANCE MATRIX FOR HOMOSCEDASTIC MLE

function _vcov!(obj::MLE{Homoscedastic}, args...)
    if getcorr(obj).method == "OIM"
        return - inv(jacobian(obj, args...))
    elseif getcorr(obj).method == "OPG"
        g = score(obj, args...)
        return g' * g
    end
end

# COVARIANCE MATRIX FOR ONE-STAGE PARAMETRIC MODELS

function _vcov!(obj::ParModel{Heteroscedastic}, args...)
    ψ     = influence(obj, args...)
    obj.V = ψ' * ψ
end

function _vcov!(obj::ParModel{<:ClusterOrCross}, args...)
    ψ     = influence(obj, args...)
    Ω     = getcorr(obj)
    obj.V = ψ' * Ω.mat * ψ
    adjfactor!(obj.V, Ω)
end

# COVARIANCE MATRIX FOR TWO-STAGE PARAMETRIC MODELS

function _vcov!(obj::TwoStageModel{Heteroscedastic}, args...)
    ψ     = influence(obj, args...)
    obj.V = ψ' * ψ
end

function _vcov!(obj::TwoStageModel{<:ClusterOrCross}, args...)
    ψ                  = influence(obj, args...)
    Ω                  = getcorr(obj)
    obj.second_stage.V = ψ' * Ω.mat * ψ
    adjfactor!(obj.second_stage.V, Ω)
end
