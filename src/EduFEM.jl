module EduFEM

export ss_beam_uniform
include("beam_solutions.jl")

export beam_mtrx
export bar_mtrx
export qe_vec
export q0_vec
include("fe_elements.jl")

export fem_solver
export stiff_mtrx
export force_member
export elem_leng
export elem_angle
export rotation_mtrx
include("fe_kernel.jl")

export postprocess
export interp_disp
export shape_func
include("fe_post.jl")

end
