# post
function postprocess(elements, nodes, u, Fₑ)

    dof = 3
    nelem = size(elements, 1)

    # nodal
    w0_ = Float64[]
    x0_ = Float64[]
    M0_ = Float64[]
    V0_ = Float64[]

    # interp
    w1_ = Float64[]
    x1_ = Float64[]

    for ielem = 1:nelem

        # element nodes
        nᵢ = elements[ielem, 1]
        nⱼ = elements[ielem, 2]

        # node coordinates
        xᵢ = nodes[nᵢ, :]
        xⱼ = nodes[nⱼ, :]

        push!(x0_, xᵢ[1], xⱼ[1])

        # index
        kᵢ = dof*(nᵢ-1)+1:dof*nᵢ
        kⱼ = dof*(nⱼ-1)+1:dof*nⱼ

        # nodal displacement
        uᵢ = u[kᵢ]
        uⱼ = u[kⱼ]

        push!(w0_, uᵢ[2], uⱼ[2])

        # interp
        xx, ww = interp_disp(xᵢ, xⱼ, uᵢ, uⱼ)

        append!(x1_, xx)
        append!(w1_, ww)

        # element force
        push!(M0_, Fₑ[ielem, 3], -Fₑ[ielem, 6])
        push!(V0_, -Fₑ[ielem, 2], Fₑ[ielem, 5])
    end

    return x0_, w0_, M0_, V0_, x1_, w1_
end

# beam displacement interpolation
function interp_disp(nodeᵢ, nodeⱼ, uᵢ, uⱼ)

    # element length
    L = elem_leng(nodeᵢ[1], nodeᵢ[2], nodeⱼ[1], nodeⱼ[2])

    dx = 1e-2
    x = collect(0.0:dx:L)

    # shape function 
    N1, N2, N3, N4 = shape_func(x, L)

    # displacement
    w = uᵢ[2] .* N1 + uᵢ[3] .* N2 + uⱼ[2] .* N3 + uⱼ[3] .* N4

    x .+= nodeᵢ[1]
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