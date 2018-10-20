#==========================================================================================#

# ONE-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj::Micromodel, corr::Homoscedastic)
end

function adjfactor!(V::Matrix, obj::Micromodel, corr::Heteroscedastic)
    if corr.adj == true
        n = nobs(obj)
        lmul!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj::Micromodel, corr::Clustered)
    if corr.adj == true
        n = corr.nc
        lmul!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj::Micromodel, corr::CrossCorrelated)
    if corr.adj == true
        n = nobs(obj)
        lmul!(n / (n - 1), V)
    end
end

#==========================================================================================#

# TWO-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj₁::Micromodel, obj₂::Micromodel, corr::Heteroscedastic)
    if corr.adj == true

        nonmissing₁ = getnonmissing(obj₁)
        nonmissing₂ = getnonmissing(obj₂)
        n           = sum(nonmissing₁ .* nonmissing₂)

        lmul!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj₁::Micromodel, obj₂::Micromodel, corr::Clustered)
    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.nonmissing .* corr₁.nonmissing .* corr₂.nonmissing
        ic    = corr.ic[touse]
        n     = length(unique(ic))

        lmul!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj₁::Micromodel, obj₂::Micromodel, corr::CrossCorrelated)
    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.nonmissing .* corr₁.nonmissing .* corr₂.nonmissing
        n     = sum(touse)

        lmul!(n / (n - 1), V)
    end
end
