
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

function _vcov(obj::MLE, corr::Homoscedastic, args...)
    if corr.method == "OIM"
        return - inv(jacobian(obj, args...))
    elseif corr.method == "OPG"
        g = score(obj, args...)
        return g' * g
    end
end

# COVARIANCE MATRIX FOR PARAMETRIC MODELS

function _vcov(obj::ParOrTwoStage, corr::Heteroscedastic, args...)
    ψ = influence(obj, args...)
    return ψ' * ψ
end

function _vcov(obj::ParOrTwoStage, corr::ClusterOrCross, args...)
    ψ = influence(obj, args...)
    return adjfactor!(ψ' * corr.mat * ψ, corr)
end
