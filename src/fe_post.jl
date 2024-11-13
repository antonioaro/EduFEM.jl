# post
function postprocess(elements, nodes, type, u, Fₑ)

    dof = 3
    nelem = length(elements)

    # nodal
    w0_ = Float64[]
    x0_ = Float64[]
    M0_ = Float64[]
    V0_ = Float64[]

    # interp
    w1_ = Float64[]
    x1_ = Float64[]

    for ielem = 1:nelem

        # element 
        n1, n2, itype, el_load = elements[ielem]
        el_type, A, I, E = type[itype]

        # node 
        x1, y1, bc1, f1 = nodes[n1]
        x2, y2, bc2, f2 = nodes[n2]

        push!(x0_, x1, x2)

        # index
        kᵢ = dof*(n1-1)+1:dof*n1
        kⱼ = dof*(n2-1)+1:dof*n2

        # nodal displacement
        uᵢ = u[kᵢ]
        uⱼ = u[kⱼ]

        push!(w0_, uᵢ[2], uⱼ[2])

        # interp
        xx, ww = interp_disp(x1, y1, x2, y2, uᵢ, uⱼ)

        append!(x1_, xx)
        append!(w1_, ww)

        # element force
        push!(M0_, -Fₑ[ielem, 3], Fₑ[ielem, 6])
        push!(V0_, -Fₑ[ielem, 2], Fₑ[ielem, 5])
    end

    return x0_, w0_, M0_, V0_, x1_, w1_
end

# beam displacement interpolation
function interp_disp(x1, y1, x2, y2, uᵢ, uⱼ)

    # element length
    L = elem_leng(x1, y1, x2, y2)

    dx = 1e-2
    x = collect(0.0:dx:L)

    # shape function 
    N1, N2, N3, N4 = shape_func(x, L)

    # displacement
    w = uᵢ[2] .* N1 + uᵢ[3] .* N2 + uⱼ[2] .* N3 + uⱼ[3] .* N4

    x .+= x1
    return x, w
end

# beam element shape function
function shape_func(x, L)
    N1 = 1 .- 3 .* (x ./ L) .^ 2 .+ 2 .* (x ./ L) .^ 3
    N2 = L .* (x ./ L .- 2 .* (x ./ L) .^ 2 .+ (x ./ L) .^ 3)
    N3 = 3 .* (x ./ L) .^ 2 .- 2 .* (x ./ L) .^ 3
    N4 = L .* (-(x ./ L) .^ 2 .+ (x ./ L) .^ 3)
    return N1, N2, N3, N4
end

# beam element shape function
function bar_shape_func(x, L)
    N1 = 1 .- (x ./ L) 
    N2 = (x ./ L )
    return N1, N2
end