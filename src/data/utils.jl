#==========================================================================================#

# COPY

copy(m::Microdata)     = Microdata([copy(getfield(m, k))     for k in fieldnames(m)]...)
deepcopy(m::Microdata) = Microdata([deepcopy(getfield(m, k)) for k in fieldnames(m)]...)

# SIZE

size(MD::Microdata)         = size(MD.mat.m)
size(MD::Microdata, dim...) = size(MD.mat.m, dim...)

#==========================================================================================#

# INTERSECTION BETWEEN A FRAME AND CORRELATION STRUCTURE

adjmsng!(msng::BitVector, corr::Homoscedastic)   = copy(corr)
adjmsng!(msng::BitVector, corr::Heteroscedastic) = copy(corr)

function adjmsng!(msng::BitVector, corr::Clustered)

    (msng == corr.msng) && (return copy(corr))

    msng   .*= corr.msng
    touse    = msng[corr.msng]
    new_ic   = corr.ic[touse]
    iter     = sort(unique(corr.ic))
    new_iter = sort(unique(new_ic))
    new_nc   = length(new_iter)
    new_mat  = corr.mat[findin(iter, new_iter), touse]

    return Clustered(corr.adj, msng, new_mat, new_ic, new_nc)
end

function adjmsng!(msng::BitVector, corr::CrossCorrelated)

    (msng == corr.msng) && (return copy(corr))

    msng  .*= corr.msng
    touse   = msng[corr.msng]
    new_mat = Symmetric(corr.mat[touse, touse])

    return CrossCorrelated(corr.adj, msng, new_mat)
end

#==========================================================================================#

# ASSIGN COLUMNS TO VARIABLE SET

function assign_columns(input::String, urterms::Terms, assign::Vector{Int})

    terms  = StatsModels.getterms(parse(input))
    output = findin(assign, findin(urterms.terms, terms))

    if any(terms .== 1) | any(terms .== +1)
        if urterms.intercept
            iszero(output) ? (output = [1]) : (output = vcat(output, 1))
        else
            throw("no intercept in model matrix for assignment")
        end
    end

    return output
end

#==========================================================================================#

# CHECK RANK

function checkrank(MD, args...)
    mat = getmatrix(MD, args)
    (rank(mat) == size(mat, 2)) || throw("model matrix does not have full rank")
end
