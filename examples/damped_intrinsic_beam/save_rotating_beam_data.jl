using Serialization: serialize

include(joinpath(@__DIR__, "elixir_rotating_beam.jl"))

output_file = isempty(ARGS) ?
              joinpath(@__DIR__, "results", "rotating_beam.jls") :
              abspath(ARGS[1])
mkpath(dirname(output_file))

node_coordinates = vec(copy(semi.cache.elements.node_coordinates))
state_shape = (12, size(semi.cache.elements.node_coordinates, 2),
               size(semi.cache.elements.node_coordinates, 3))
states = [Array(reshape(state, state_shape)) for state in sol.u]
steady_states = [collect(steady_rotating_solution(x)) for x in node_coordinates]

data = (times = collect(sol.t),
        node_coordinates,
        states,
        flexibility_matrix = Matrix(flexibility_matrix),
        initial_curvature = collect(equations_hyperbolic.initial_curvature),
        energy_history,
        terminal_angular_speed,
        ramp_duration,
        steady_states)
serialize(output_file, data)
println("wrote ", output_file)
