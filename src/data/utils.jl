#==========================================================================================#

# COPY

copy(m::Microdata)     = Microdata([copy(getfield(m, k))     for k in fieldnames(m)]...)
deepcopy(m::Microdata) = Microdata([deepcopy(getfield(m, k)) for k in fieldnames(m)]...)

#==========================================================================================#

# INTERSECTION BETWEEN A FRAME AND CORRELATION STRUCTURE

adjmissing!(nonmissing::BitVector, corr::Homoscedastic)   = copy(corr)
adjmissing!(nonmissing::BitVector, corr::Heteroscedastic) = copy(corr)

function adjmissing!(nonmissing::BitVector, corr::Clustered)

    (nonmissing == corr.nonmissing) && (return copy(corr))

    nonmissing .*= corr.nonmissing
    touse        = nonmissing[corr.nonmissing]
    new_ic       = corr.ic[touse]
    iter         = sort(unique(corr.ic))
    new_iter     = sort(unique(new_ic))
    new_nc       = length(new_iter)
    new_mat      = corr.mat[findall((in)(new_iter), iter), touse]

    return Clustered(corr.adj, nonmissing, new_mat, new_ic, new_nc)
end

function adjmissing!(nonmissing::BitVector, corr::CrossCorrelated)

    (nonmissing == corr.nonmissing) && (return copy(corr))

    nonmissing .*= corr.nonmissing
    touse        = nonmissing[corr.nonmissing]
    new_mat      = Symmetric(corr.mat[touse, touse])

    return CrossCorrelated(corr.adj, nonmissing, new_mat)
end

#==========================================================================================#

# ASSIGN COLUMNS TO VARIABLE SET

function assign_columns(input::String, urterms::StatsModels.Terms, assign::Vector{Int})

    formula = @eval @formula $nothing ~ $(Meta.parse(input))
    terms   = StatsModels.Terms(formula)
    output  = findall((in)(findall((in)(terms.terms), urterms.terms)), assign)

    if occursin("1", input)
        iszero(output) ? (output = [1]) : (output = vcat(output, 1))
    end

    return output
end

#==========================================================================================#

# CHECK RANK

function checkrank(MD, args...)
    mat = getmatrix(MD, args)
    (rank(mat) == size(mat, 2)) || throw("model matrix does not have full rank")
end
