#==========================================================================================#

# INTERSECTION BETWEEN A FRAME AND CORRELATION STRUCTURE

function _adjmsng!(msng::BitVector, corr::Homoscedastic, makecopy::Bool)
    return (makecopy ? copy(corr) : corr)
end

function _adjmsng!(msng::BitVector, corr::Heteroscedastic, makecopy::Bool)
    return (makecopy ? copy(corr) : corr)
end

function _adjmsng!(msng::BitVector, corr::Clustered, makecopy::Bool)

    (msng == corr.msng) && (return (makecopy ? copy(corr) : corr))

    intersection = msng .* corr.msng
    msng[:]      = intersection
    touse        = intersection[corr.msng]
    new_ic       = (makecopy ? copy(corr.ic[touse]) : view(corr.ic, touse))
    new_mat      = (makecopy ? copy(corr.mat[touse, touse]) : view(corr.mat, touse, touse))

    return Clustered(intersection, new_mat, new_ic, length(unique(new_ic)))
end

function _adjmsng!(msng::BitVector, corr::CrossCorrelated, makecopy::Bool)

    (msng == corr.msng) && (return (makecopy ? copy(corr) : corr))

    intersection = msng .* corr.msng
    msng[:]      = intersection
    touse        = intersection[corr.msng]
    new_mat      = (makecopy ? copy(corr.mat[touse, touse]) : view(corr.mat, touse, touse))

    return CrossCorrelated(intersection, new_mat)
end

#==========================================================================================#

# ASSIGN COLUMNS TO VARIABLE SET

function assign_columns(input::String, urterms::Terms, assign::Vector{Int})

    terms  = DataFrames.getterms(parse(input))
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
