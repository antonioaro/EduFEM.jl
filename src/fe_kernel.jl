# solver
function fem_solver(elements, nodes, type, boundary_condition, nodal_force, element_force)

    dof = 3
    nnode = size(nodes, 1)
    ndof = nnode * dof

    # displacement and force vector
    u = zeros(ndof)

    # boundary condition
    dof_d = get_dof_d(nodes, boundary_condition)
    dof_f = filter(x -> x in 1:ndof && !(x in dof_d), 1:ndof)

    # ecuación de equilibrio
    K, F, Q = get_system(ndof, elements, nodes, type, nodal_force, element_force)

    # solution
    u[dof_f] = K[dof_f, dof_f] \ (F[dof_f] + Q[dof_f])

    # reaction force
    F[dof_d] = K[dof_d, dof_d] * u[dof_d] + K[dof_d, dof_f] * u[dof_f] - Q[dof_d]

    # element force member
    Fₑ = force_member(elements, nodes, type, element_force, u)

    return u, F, Fₑ
end

# global stiffness matrix
function get_system(ndof, elements, nodes, type, nodal_force, element_force)

    dof = 3
    nelem = length(elements)

    K = zeros(ndof, ndof)
    F = zeros(ndof, 1)
    Q = zeros(ndof, 1)

    for ielem = 1:nelem

        # element 
        n1, n2, itype, el_load = elements[ielem]
        el_type, A, I, E = type[itype]

        # node 
        x1, y1, bc1, f1 = nodes[n1]
        x2, y2, bc2, f2 = nodes[n2]

        # element length
        leng = elem_leng(x1, y1, x2, y2)

        # element angle
        αₑ = elem_angle(x1, y1, x2, y2)
        Ld = rotation_mtrx(αₑ)

        # index
        kᵢ = dof*(n1-1)+1:dof*n1
        kⱼ = dof*(n2-1)+1:dof*n2

        # element matrix
        if el_type == "bar"
            Ke = bar_mtrx(E, A, leng)
        elseif el_type == "beam"
            Ke = beam_mtrx(E, A, I, leng)
        else
            error("Not valid element")
        end

        # assemble
        index = vcat(kᵢ, kⱼ)
        K[index, index] += Ld * Ke * Ld'

        # element force vector        
        if el_load !== 0
            label, qx, qz, e0 = element_force[el_load]

            Qₑ = qe_vec(qx, qz, leng)
            Q₀ = q0_vec(e0, E, A)
            Q[index] += Ld * (Qₑ + Q₀)
        end
    end

    # nodal load
    nnode = length(nodes)
    for inode = 1:nnode
        fnode = Int(nodes[inode][4])
        if fnode !== 0
            label, fx, fy, Mz = nodal_force[fnode]
            idx = dof*(inode-1)+1:dof*inode
            F[idx] .= [fx, fy, Mz]
        end
    end

    return K, F, Q
end

# element force member
function force_member(elements, nodes, type, element_force, u)

    dof = 3
    edof = 6
    nelem = length(elements)
    Fₑ = zeros(nelem, edof)

    for ielem = 1:nelem

        # element 
        n1, n2, itype, el_load = elements[ielem]
        el_type, A, I, E = type[itype]

        # node 
        x1, y1, bc1, f1 = nodes[n1]
        x2, y2, bc2, f2 = nodes[n2]

        # element length
        leng = elem_leng(x1, y1, x2, y2)

        # element angle
        αₑ = elem_angle(x1, y1, x2, y2)
        Ld = rotation_mtrx(αₑ)

        # index
        kᵢ = dof*(n1-1)+1:dof*n1
        kⱼ = dof*(n2-1)+1:dof*n2

        # element matrix
        if el_type == "bar"
            Ke = bar_mtrx(E, A, leng)
        elseif el_type == "beam"
            Ke = beam_mtrx(E, A, I, leng)
        else
            error("Not valid element")
        end

        # element displacement
        index = vcat(kᵢ, kⱼ)
        uₑ = Ld' * u[index]

        # element force
        Fₑ[ielem, :] = Ke * uₑ

        # element force vector    
        # if el_load !== 0
        #     label, qx, qz, e0 = element_force[el_load]

        #     Q₀ = Q0_vec(e0, E, A)
        #     Fₑ[ielem, :] -= Q₀
        # end
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

# get dof with bc
function get_dof_d(nodes, boundary_condition)
    dof = 3
    dof_d = Vector{Int}()
    nnodes = length(nodes)

    for inode = 1:nnodes
        ibc = Int(nodes[inode][3])

        if ibc != 0
            bc = boundary_condition[ibc]

            for idof = 1:dof
                if !isnan(bc[idof+1])
                    index = dof * (inode - 1) + idof
                    push!(dof_d, index)
                end
            end
        end
    end

    return dof_d
end