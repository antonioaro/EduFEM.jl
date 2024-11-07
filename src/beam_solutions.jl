# simply supported beam under uniform load
function ss_beam_uniform(E, I, L, x, p)
    y = p .* x ./ (24 * E * I) .* (x .^ 3 - 2 * L .* x .^ 2 .+ L^3)
    M = p / 2 .* x .* (L .- x)
    V = p .* (L / 2 .- x)
    return y, M, V
end