module EduFEM

export ss_beam_uniform
include("beam_solutions.jl")

export beam_mtrx
export bar_mtrx
export qe_vec
export q0_vec
export bar_shape_func
include("fe_elements.jl")

export fem_solver
export get_system
export force_member
export elem_leng
export elem_angle
export rotation_mtrx
export get_dof_d
include("fe_kernel.jl")

export postprocess
export interp_disp
export shape_func
include("fe_post.jl")

export fe_mesh
export equivalence
include("fe_mesh.jl")

end
