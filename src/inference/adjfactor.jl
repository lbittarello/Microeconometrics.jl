#==========================================================================================#

# ONE-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(V::Matrix, obj::Micromodel{Heteroscedastic})

    Ω = getcorr(obj)

    if Ω.adj == true
        n = nobs(obj)
        scale!(n / (n - 1), V)
    end
end


function adjfactor!(V::Matrix, obj::Micromodel{CrossCorrelated})

    Ω = getcorr(obj)

    if Ω.adj == true
        n = nobs(obj)
        scale!(n / (n - 1), V)
    end
end

function adjfactor!(V::Matrix, obj::Micromodel{Clustered})

    Ω = getcorr(obj)

    if Ω.adj == true
        n = getcorr(obj).nc
        scale!(n / (n - 1), V)
    end
end

#==========================================================================================#

# TWO-SAMPLE FINITE-SAMPLE ADJUSTMENT

function adjfactor!(
        V::Matrix, obj₁::Micromodel{T}, obj₂::Micromodel{T}, corr::T
    ) where {T <: Heteroscedastic}

    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.msng .* corr₁.msng .* corr₂.msng
        n     = sum(touse)

        scale!(n / (n - 1), V)
    end
end

function adjfactor!(
        V::Matrix, obj₁::Micromodel{T}, obj₂::Micromodel{T}, corr::T
    ) where {T <: CrossCorrelated}

    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.msng .* corr₁.msng .* corr₂.msng
        n     = sum(touse)

        scale!(n / (n - 1), V)
    end
end

function adjfactor!(
        V::Matrix, obj₁::Micromodel{T}, obj₂::Micromodel{T}, corr::T
    ) where {T <: Clustered}

    if corr.adj == true

        corr₁ = getcorr(obj₁)
        corr₂ = getcorr(obj₂)
        touse = corr.msng .* corr₁.msng .* corr₂.msng
        ic    = corr.ic[touse]
        n     = length(unique(ic))

        scale!(n / (n - 1), V)
    end
end
