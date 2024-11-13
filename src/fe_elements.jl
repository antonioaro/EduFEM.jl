# beam element
function beam_mtrx(E, A, I, L)
    # Matriz de rigidez para una viga en 2D
    Ke = [
        E*A/L 0 0 -E*A/L 0 0;
        0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
        0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
        -E*A/L 0 0 E*A/L 0 0;
        0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
        0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L
    ]

    return Ke
end

# bar element
function bar_mtrx(E, A, L)
    Ke = E * A / L * [1 0 0 -1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        -1 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0]
    return Ke
end

# element load
function qe_vec(Qₓ, Qz, L)
    Qe = [Qₓ * L / 2;
        Qz * L / 2;
        Qz * L^2 / 12;
        Qₓ * L / 2;
        Qz * L / 2;
        -Qz * L^2 / 12]
    return Qe
end

# element load initial strain
function q0_vec(ε₀, E, A)
    Q0 = E * A * ε₀ .* [-1.0; 0; 0; 1.0; 0; 0]
    return Q0
end

# beam element shape function
function bar_shape_func(x, L)
    N1 = 1 .- (x ./ L)
    N2 = (x ./ L)
    return N1, N2
end


