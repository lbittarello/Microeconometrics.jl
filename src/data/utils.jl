#==========================================================================================#

# INTERSECTION BETWEEN A FRAME AND CORRELATION STRUCTURE

adjmsng!(msng::BitVector, corr::Homoscedastic)   = copy(corr)
adjmsng!(msng::BitVector, corr::Heteroscedastic) = copy(corr)

function adjmsng!(msng::BitVector, corr::Clustered)

    (msng == corr.msng) && (return copy(corr))

    intersect = msng .* corr.msng
    msng     .= intersect
    touse     = intersect[corr.msng]
    new_ic    = copy(corr.ic[touse])
    new_nc    = length(unique(new_ic))
    new_mat   = copy(corr.mat[touse, touse])

    return Clustered(corr.adj, intersect, new_mat, new_ic, new_nc)
end

function adjmsng!(msng::BitVector, corr::CrossCorrelated)

    (msng == corr.msng) && (return copy(corr))

    intersect = msng .* corr.msng
    msng     .= intersect
    touse     = intersect[corr.msng]
    new_mat   = copy(corr.mat[touse, touse])

    return CrossCorrelated(corr.adj, intersect, new_mat)
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
