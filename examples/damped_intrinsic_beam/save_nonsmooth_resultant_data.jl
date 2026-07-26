using Serialization: serialize

include(joinpath(@__DIR__, "elixir_nonsmooth_resultant.jl"))

output_file = isempty(ARGS) ?
              joinpath(@__DIR__, "results", "nonsmooth_resultant.jls") :
              abspath(ARGS[1])
mkpath(dirname(output_file))

node_coordinates = vec(copy(semi.cache.elements.node_coordinates))
state_shape = (12, size(semi.cache.elements.node_coordinates, 2),
               size(semi.cache.elements.node_coordinates, 3))
states = [Array(reshape(state, state_shape)) for state in sol.u]

data = (times = collect(sol.t),
        node_coordinates,
        states,
        flexibility_matrix = Matrix(flexibility_matrix),
        initial_curvature = collect(equations_hyperbolic.initial_curvature),
        energy_history,
        maximum_axial_resultant,
        maximum_axial_strain)
serialize(output_file, data)
println("wrote ", output_file)
