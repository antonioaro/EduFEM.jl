# solver
function fem_solver(elements, nodes, type, E, A, I, Qx, Qz, ε₀, dof_d)

    dof = 3
    nnode = size(nodes, 1)
    ndof = nnode * dof

    # displacement and force vector
    u = zeros(ndof)
    f = zeros(ndof)

    # boundary condition
    dof_f = filter(x -> x in 1:ndof && !(x in dof_d), 1:ndof)

    # ecuación de equilibrio
    K, Q = stiff_mtrx(ndof, elements, nodes, type, E, A, I, Qx, Qz, ε₀)

    # solution
    u[dof_f] = K[dof_f, dof_f] \ (f[dof_f] + Q[dof_f])

    # reaction force
    f[dof_d] = K[dof_d, dof_d] * u[dof_d] + K[dof_d, dof_f] * u[dof_f] - Q[dof_d]

    # element force member
    Fₑ = force_member(elements, nodes, type, E, A, I, ε₀, u)

    return u, f, Fₑ
end

# global stiffness matrix
function stiff_mtrx(ndof, elements, nodes, type, E, A, I, Qx, Qz, ε₀)

    dof=3

    K = zeros(ndof, ndof)
    Q = zeros(ndof, 1)

    nelem = size(elements, 1)

    for ielem = 1:nelem

        # element nodes
        nᵢ = elements[ielem, 1]
        nⱼ = elements[ielem, 2]

        # node coordinates
        xᵢ = nodes[nᵢ, :]
        xⱼ = nodes[nⱼ, :]

        # index
        kᵢ = dof*(nᵢ-1)+1:dof*nᵢ
        kⱼ = dof*(nⱼ-1)+1:dof*nⱼ

        # element length
        Lₑ = elem_leng(xᵢ[1], xᵢ[2], xⱼ[1], xⱼ[2])

        # element angle
        αₑ = elem_angle(xᵢ[1], xᵢ[2], xⱼ[1], xⱼ[2])
        Ld = rotation_mtrx(αₑ)

        # element matrix
        if type[ielem] == "bar"
            Kₑ = bar_mtrx(E[ielem], A[ielem], Lₑ)
        elseif type[ielem] == "beam"
            Kₑ = beam_mtrx(E[ielem], A[ielem], I[ielem], Lₑ)
        else
            error("Not valid element")
        end

        # assemble
        index = vcat(kᵢ, kⱼ)
        K[index, index] += Ld * Kₑ * Ld'

        # element force vector
        Qₑ = qe_vec(Qx[ielem], Qz[ielem], Lₑ)
        Q₀ = q0_vec(ε₀[ielem], E[ielem], A[ielem])
        Q[index] += Ld * (Qₑ + Q₀)
    end

    return K, Q
end

# element force member
function force_member(elements, nodes, type, E, A, I, ε₀, u)

    dof = 3
    edof = 6
    nelem = size(elements, 1)
    Fₑ = zeros(nelem, edof)

    for ielem = 1:nelem
        # element nodes
        nᵢ = elements[ielem, 1]
        nⱼ = elements[ielem, 2]

        # node coordinates
        xᵢ = nodes[nᵢ, :]
        xⱼ = nodes[nⱼ, :]

        # index
        kᵢ = dof*(nᵢ-1)+1:dof*nᵢ
        kⱼ = dof*(nⱼ-1)+1:dof*nⱼ

        # element length
        Lₑ = elem_leng(xᵢ[1], xᵢ[2], xⱼ[1], xⱼ[2])

        # element angle
        αₑ = elem_angle(xᵢ[1], xᵢ[2], xⱼ[1], xⱼ[2])
        Ld = rotation_mtrx(αₑ)

        # element matrix
        if type[ielem] == "bar"
            Kₑ = bar_mtrx(E[ielem], A[ielem], Lₑ)
        elseif type[ielem] == "beam"
            Kₑ = beam_mtrx(E[ielem], A[ielem], I[ielem], Lₑ)
        else
            error("Not valid element")
        end

        # element displacement
        index = vcat(kᵢ, kⱼ)
        uₑ = Ld' * u[index]

        # element force
        Fₑ[ielem, :] = Kₑ * uₑ

        if ε₀[ielem] != 0
            Q₀ = Q0_vec(ε₀[ielem], E[ielem], A[ielem])
            Fₑ[ielem, :] -= Q₀
        end
    end
    return Fₑ
end

# element length
function elem_leng(x1, y1, x2, y2)
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    return L
end

# element angle
function elem_angle(x1, y1, x2, y2)
    αₑ = atan((y2 - y1) / (x2 - x1))
    return αₑ
end

# rotation matrix
function rotation_mtrx(a)
    Ld = [cos(a) -sin(a) 0 0 0 0;
        sin(a) cos(a) 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 cos(a) -sin(a) 0;
        0 0 0 sin(a) cos(a) 0;
        0 0 0 0 0 1]
    return Ld
end