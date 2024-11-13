function fe_geom(elements, nodes)

    nelem = size(elements, 1)

    # nodal
    x_ = Float64[]

    for ielem = 1:nelem

        # element nodes
        nᵢ = elements[ielem, 1]
        nⱼ = elements[ielem, 2]

        # node coordinates
        xᵢ = nodes[nᵢ, :]
        xⱼ = nodes[nⱼ, :]

        push!(x_, xᵢ[1], xⱼ[1])        
    end

    return x_
end


function fe_mesh(elements, nodes, Ne=0, Le=0.0)

    elements_ = Vector{Vector{Int}}()
    nodes_ = Vector{Vector{Float64}}()

    nelem = size(elements, 1)

    for ielem = 1:nelem

        # element 
        n1, n2, type, el_load = elements[ielem]

        # node 
        x1, y1, bc1, f1 = nodes[n1]
        x2, y2, bc2, f2 = nodes[n2]

        # element length
        leng = elem_leng(x1, y1, x2, y2)

        if Le != 0.0
            Ne = Int(floor(leng / Le))
        end

        # shape function
        x = range(0.0, leng, length=Ne + 1)
        N1, N2 = bar_shape_func(x, leng)

        # new_nodes
        offset = size(nodes_, 1)
        new_nodes = [[N1[i] * x1 + N2[i] * x2, N1[i] * y1 + N2[i] * y2, 0, 0] for i in 1:length(x)]

        # boundary conditions
        new_nodes[1][3] = bc1
        new_nodes[1][4] = f1
        new_nodes[end][3] = bc2
        new_nodes[end][4] = f2
        append!(nodes_, new_nodes)

        # new elements
        new_elements = [[offset + i, offset + i + 1, type, el_load] for i in 1:Ne]
        append!(elements_, new_elements)
    end

    elements_, nodes_ = equivalence(elements_, nodes_)

    return elements_, nodes_
end

function equivalence(elements, nodes)
    nnode = length(nodes)
    nelem = length(elements)

    duplicated_node = []
    by_node = []

    for inode = 1:nnode

        if inode in duplicated_node
            continue
        end

        x1, y1, bc1, f1 = nodes[inode]

        for jnode = 1:nnode
            if inode != jnode
                x2, y2, bc2, f2 = nodes[jnode]

                # check for duplicated node
                if isapprox(x1, x2, atol=1e-6) && isapprox(y1, y2, atol=1e-6)
                    push!(duplicated_node, jnode)
                    push!(by_node, inode)
                end
            end
        end
    end

    # not duplicated nodes
    unique_node = setdiff(1:nnode, duplicated_node)
    nodes_ = [nodes[i] for i in unique_node]

    # replace node in element   
    for ielem = 1:nelem
        n1, n2, type, el_load = elements[ielem]

        # Reemplazar nodos duplicados en los elementos
        if n1 in duplicated_node
            i1 = findfirst(x -> x == n1, duplicated_node)
            n1 = by_node[i1]
        end

        if n2 in duplicated_node
            i2 = findfirst(x -> x == n2, duplicated_node)
            n2 = by_node[i2]  # Reemplazar por el nodo original
        end

        # Actualizar los elementos con los nuevos índices
        elements[ielem][1] = findfirst(x -> x == n1, unique_node)
        elements[ielem][2] = findfirst(x -> x == n2, unique_node)
    end

    return elements, nodes_
end