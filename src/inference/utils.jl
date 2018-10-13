#==========================================================================================#

# COPY

copy(corr::Homoscedastic)   = Homoscedastic(corr.adj, corr.method)
copy(corr::Heteroscedastic) = Heteroscedastic(corr.adj)

function copy(corr::CrossCorrelated)
    CrossCorrelated(corr.adj, copy(corr.nonmissing), copy(corr.mat))
end

function copy(corr::Clustered)
    Clustered(corr.adj, copy(corr.nonmissing), copy(corr.mat), copy(corr.ic), copy(corr.nc))
end

# EQUALITY

function Base.isequal(corr₁::C, corr₂::C) where {C <: CorrStructure}

    output = true

    for k in fieldnames(typeof(corr₁))
        output *= isequal(getfield(corr₁, k), getfield(corr₂, k))
    end

    return output
end

Base.isequal(corr₁::CorrStructure, corr₂::CorrStructure) = false

#==========================================================================================#

# DISTANCE ON THE GLOBE

havd(x::Float64) = 0.5 * (1.0 - cosd(x))

function geodistance(y0::Float64, x0::Float64, y1::Float64, x1::Float64)
    h = havd(y1 - y0) + havd(x1 - x0) * cosd(y0) * cosd(y1)
    return 12742.0 * asin(min(1.0, sqrt(h)))
end

#==========================================================================================#

# KERNELS

bartlett(x::Real)        = bartlett(float(x))
bartlett(x::Float64)     = (abs(x) >= 1.0) ? 0.0 : (1.0 - abs(x))
truncated(x::Real)       = truncated(float(x))
truncated(x::Float64)    = (abs(x) >= 1.0) ? 0.0 : 1.0
tukeyhanning(x::Real)    = tukeyhanning(float(x))
tukeyhanning(x::Float64) = (abs(x) >= 1.0) ? 0.0 : 0.5 * (1.0 + cospi(x))

parzen(x::Real) = parzen(float(x))

function parzen(x::Float64)
    if abs(x) <= 0.5
        return 1.0 - 6.0 * abs2(x) + 6 * abs(x)^3
    elseif abs(x) <= 1.0
        return 2.0 * (1.0 - abs(x))^3
    else
        return 0.0
    end
end

gallant(x::Real)    = parzen(float(x))
gallant(x::Float64) = parzen(x)
