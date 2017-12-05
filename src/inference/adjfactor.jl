#==========================================================================================#

# ONE-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj::Micromodel, corr::Heteroscedastic)
    if corr.adj == true
        n = nobs(obj)
        scale!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj::Micromodel, corr::Clustered)
    if corr.adj == true
        n = corr.nc
        scale!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj::Micromodel, corr::CrossCorrelated)
    if corr.adj == true
        n = nobs(obj)
        scale!(n / (n - 1), V)
    end
end

#==========================================================================================#

# TWO-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj₁::Micromodel, obj₂::Micromodel, corr::Heteroscedastic)
    if corr.adj == true

        msng₁ = getmsng(obj₁)
        msng₂ = getmsng(obj₂)
        n     = sum(msng₁ .* msng₂)
        
        scale!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj₁::Micromodel, obj₂::Micromodel, corr::Clustered)
    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.msng .* corr₁.msng .* corr₂.msng
        ic    = corr.ic[touse]
        n     = length(unique(ic))

        scale!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj₁::Micromodel, obj₂::Micromodel, corr::CrossCorrelated)
    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.msng .* corr₁.msng .* corr₂.msng
        n     = sum(touse)

        scale!(n / (n - 1), V)
    end
end
