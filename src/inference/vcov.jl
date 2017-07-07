#==========================================================================================#

# CONTRIBUTION OF FIRST-STAGE ESTIMATION TO SECOND-STAGE INFLUENCE FUNCTION

function crossinfluence(obj::TwoStageModel, obj₁::ParModel, args...)
    return influence(obj₁, args...) * crossjacobian(obj, args...)'
end

#==========================================================================================#

# EXPECTED OUTER PRODUCT OF THE SCORE VECTOR × NUMBER OF OBSERVATIONS
# INVERSE OF OPTIMAL WEIGHT MATRIX FOR GMM

function _opg(obj::ParModel, corr::Heteroscedastic, args...)
    ψ = score(obj, args...)
    return ψ' * ψ
end

function _opg(obj::ParModel, corr::ClusterOrCross, args...)
    ψ = score(obj, args...)
    return ψ' * corr.mat * ψ
end

function _opg(obj::TwoStageModel, corr::Heteroscedastic, args...)
    ψ = score(obj, args...) + crossinfluence(obj, first_stage(obj), args...)
    return ψ' * ψ
end

function _opg(obj::TwoStageModel, corr::ClusterOrCross, args...)
    ψ = score(obj, args...) + crossinfluence(obj, first_stage(obj), args...)
    return ψ' * corr.mat * ψ
end

#==========================================================================================#

# INFLUENCE FUNCTION

function influence(obj::ParModel, args...)
    s = score(obj, args...)
    j = jacobian(obj, args...)
    if obj.method == "Unadjusted"
        return scale!(- 1.0, (s * j) / (j' * j))
    else
        return scale!(- 1.0, s / j')
    end
end

function influence(obj::TwoStageModel, obj₁::Micromodel, obj₂::ParModel, args...)
    s = score(obj, args...) + crossinfluence(obj, obj₁, args...)
    j = jacobian(obj, args...)
    if obj₂.method == "Unadjusted"
        return scale!(- 1.0, (s * j) / (j' * j))
    else
        return scale!(- 1.0, s / j')
    end
end

function influence(obj::TwoStageModel, args...)
    return influence(obj, first_stage(obj), second_stage(obj), args...)
end

#==========================================================================================#

# COVARIANCE MATRIX FOR HOMOSCEDASTIC MLE

function _vcov(obj::ParModel, corr::Homoscedastic, args...)
    if obj.method == "MLE"
        if corr.method == "OIM"
            return - inv(jacobian(obj, args...))
        elseif corr.method == "OPG"
            g = score(obj, args...)
            return g' * g
        end
    else
        throw("homoscecastic errors are only available for OLS and MLE")
    end
end

# COVARIANCE MATRIX FOR PARAMETRIC MODELS

function _vcov(obj::ParOrTwoStage, corr::Heteroscedastic, args...)
    ψ = influence(obj, args...)
    return ψ' * ψ
end

function _vcov(obj::ParOrTwoStage, corr::ClusterOrCross, args...)
    ψ = influence(obj, args...)
    return _adjcluster!(ψ' * corr.mat * ψ, corr)
end
