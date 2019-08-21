#==========================================================================================#

# ONE-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj::ParModel, corr::Homoscedastic)
end

function adjfactor!(V::Matrix, obj::ParModel, corr::Heteroscedastic)
    (corr.adj == true) && (n = nobs(obj) ; lmul!(n / (n - 1), V))
end

function adjfactor!(V::Matrix, obj::ParModel, corr::Clustered)
    (corr.adj == true) && (n = corr.nc ; lmul!(n / (n - 1), V))
end

function adjfactor!(V::Matrix, obj::ParModel, corr::CrossCorrelated)
    (corr.adj == true) && (n = nobs(obj) ; lmul!(n / (n - 1), V))
end

#==========================================================================================#

# TWO-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj₁::ParModel, obj₂::ParModel, corr::Heteroscedastic)
    if corr.adj == true

        nonmissing₁ = getnonmissing(obj₁)
        nonmissing₂ = getnonmissing(obj₂)
        n           = sum(nonmissing₁ .* nonmissing₂)

        lmul!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj₁::ParModel, obj₂::ParModel, corr::Clustered)
    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.nonmissing .* corr₁.nonmissing .* corr₂.nonmissing
        ic    = corr.ic[touse]
        n     = length(unique(ic))

        lmul!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj₁::ParModel, obj₂::ParModel, corr::CrossCorrelated)
    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.nonmissing .* corr₁.nonmissing .* corr₂.nonmissing
        n     = sum(touse)

        lmul!(n / (n - 1), V)
    end
end
