using Dates
using LinearAlgebra: eigvals
using Printf

# Load the benchmark once. The short setup solve keeps the elixir directly
# executable while avoiding an additional long fixed-frequency calculation
# before the continuation sweep.
ENV["CANTILEVER_SETTLING_CYCLES"] = "1"
ENV["CANTILEVER_MEASUREMENT_CYCLES"] = "2"
ENV["CANTILEVER_RAMP_CYCLES"] = "0"
ENV["CANTILEVER_SAMPLES_PER_CYCLE"] =
    get(ENV, "CANTILEVER_SAMPLES_PER_CYCLE", "96")
ENV["CANTILEVER_POLYDEG"] = get(ENV, "CANTILEVER_POLYDEG", "4")
ENV["CANTILEVER_REFINEMENT_LEVEL"] =
    get(ENV, "CANTILEVER_REFINEMENT_LEVEL", "1")
Base.include(@__MODULE__,
             joinpath(@__DIR__, "elixir_base_excited_cantilever.jl"))

function dense_finite_difference_jacobian(state)
    dimension = length(state)
    jacobian = zeros(dimension, dimension)
    base_derivative = similar(state)
    perturbed_derivative = similar(state)
    perturbed_state = copy(state)
    cantilever_rhs!(base_derivative, state, rhs_parameters, 0.0)
    for column in eachindex(state)
        step = sqrt(eps(Float64)) * max(1.0, abs(state[column]))
        perturbed_state[column] += step
        cantilever_rhs!(perturbed_derivative, perturbed_state,
                        rhs_parameters, 0.0)
        perturbed_state[column] = state[column]
        @views jacobian[:, column] .=
            (perturbed_derivative .- base_derivative) ./ step
    end
    return jacobian
end

function first_discrete_mode(state)
    eigenvalues = eigvals(dense_finite_difference_jacobian(state))
    oscillatory = sort([value for value in eigenvalues if imag(value) > 0.1];
                       by = imag)
    isempty(oscillatory) &&
        error("no oscillatory eigenvalue was found")
    return first(oscillatory)
end

function response_metrics(states)
    tips = [reconstruct_tip(state) for state in states]
    samples = samples_per_cycle
    last_indices = (length(tips) - samples):length(tips)
    previous_indices = ((length(tips) - 2 * samples):
                        (length(tips) - samples))
    last_tip = tips[last_indices]
    previous_tip = tips[previous_indices]
    return (;
            transverse_peak = maximum(abs(value[2]) for value in last_tip),
            longitudinal_minimum = minimum(value[1] for value in last_tip),
            rotation_peak = maximum(abs(value[3]) for value in last_tip),
            periodicity_error = maximum(maximum(abs,
                                                 last_tip[index] -
                                                 previous_tip[index])
                                        for index in eachindex(last_tip)))
end

function solve_cycles(initial_state, cycles, frequency;
                      ramp_cycles = 0.0)
    root_velocity.omega = 2.0 * pi * frequency * FAROKHI_T
    local_period = 2.0 * pi / root_velocity.omega
    root_velocity.ramp_duration = ramp_cycles * local_period
    measurement_cycles = 2
    measurement_start = (cycles - measurement_cycles) * local_period
    save_times = range(measurement_start, cycles * local_period;
                       length = measurement_cycles *
                                samples_per_cycle + 1)
    local_problem = remake(problem;
                           u0 = copy(initial_state),
                           tspan = (0.0, cycles * local_period))
    local_solution = solve(local_problem, algorithm;
                           reltol = parse(Float64,
                                          get(ENV,
                                              "CANTILEVER_SWEEP_RELTOL",
                                              "1e-4")),
                           abstol = parse(Float64,
                                          get(ENV,
                                              "CANTILEVER_SWEEP_ABSTOL",
                                              "1e-6")),
                           dtmax = local_period /
                                   parse(Float64,
                                         get(ENV,
                                             "CANTILEVER_STEPS_PER_PERIOD",
                                             "32")),
                           saveat = save_times,
                           save_start = false,
                           save_everystep = false,
                           maxiters = 10^7)
    @assert SciMLBase.successful_retcode(local_solution)
    return local_solution, response_metrics(local_solution.u)
end

normalized_frequencies = parse.(Float64,
                                split(get(ENV,
                                          "CANTILEVER_SWEEP_NORMALIZED_FREQUENCIES",
                                          "0.9003,0.9517,0.9774,0.9928,1.0031,1.0083,1.0124,1.0165,1.0206"),
                                      ','))
first_cycles = parse(Int, get(ENV, "CANTILEVER_SWEEP_FIRST_CYCLES", "40"))
continuation_cycles = parse(Int,
                            get(ENV, "CANTILEVER_SWEEP_CONTINUATION_CYCLES",
                                "18"))
additional_cycles = parse(Int,
                          get(ENV, "CANTILEVER_SWEEP_ADDITIONAL_CYCLES",
                              "10"))
maximum_cycles = parse(Int,
                       get(ENV, "CANTILEVER_SWEEP_MAXIMUM_CYCLES", "60"))
periodicity_tolerance = parse(Float64,
                              get(ENV,
                                  "CANTILEVER_SWEEP_PERIODICITY_TOLERANCE",
                                  "2e-3"))

mode = first_discrete_mode(split_problem.u0)
linear_frequency_hz = imag(mode) / (2.0 * pi * FAROKHI_T)
linear_damping_ratio = -real(mode) / abs(mode)
frequencies = normalized_frequencies .* linear_frequency_hz
@printf("Discrete first mode: lambda = %.8e %+.8ei, f1 = %.8f Hz, zeta = %.6e\n",
        real(mode), imag(mode), linear_frequency_hz, linear_damping_ratio)

rows = NamedTuple[]
let continuation_state = copy(split_problem.u0)
    for (index, frequency) in enumerate(frequencies)
        cycles = index == 1 ? first_cycles : continuation_cycles
        ramp = index == 1 ?
               parse(Float64,
                     get(ENV, "CANTILEVER_SWEEP_RAMP_CYCLES", "8")) :
               0.0
        local_solution, metrics = solve_cycles(continuation_state, cycles,
                                               frequency;
                                               ramp_cycles = ramp)
        continuation_state = copy(local_solution.u[end])
        total_cycles = cycles

        while metrics.periodicity_error > periodicity_tolerance &&
              total_cycles < maximum_cycles
            extra_solution, metrics = solve_cycles(continuation_state,
                                                    additional_cycles,
                                                    frequency)
            continuation_state = copy(extra_solution.u[end])
            total_cycles += additional_cycles
            local_solution = extra_solution
        end

        row = (;
               acceleration_rms_g,
               frequency_hz = frequency,
               frequency_over_discrete_f1 = frequency / linear_frequency_hz,
               omega = root_velocity.omega,
               transverse_peak = metrics.transverse_peak,
               longitudinal_minimum = metrics.longitudinal_minimum,
               rotation_peak = metrics.rotation_peak,
               periodicity_error = metrics.periodicity_error,
               cycles = total_cycles,
               accepted_steps = local_solution.destats.naccept,
               rejected_steps = local_solution.destats.nreject)
        push!(rows, row)
        @printf("f=%7.4f Hz  f/f1=%8.5f  |w|=%9.6f  u_min=%9.6f  |psi|=%9.6f  periodic=%8.2e  cycles=%d\n",
                row.frequency_hz, row.frequency_over_discrete_f1,
                row.transverse_peak, row.longitudinal_minimum,
                row.rotation_peak, row.periodicity_error, row.cycles)
    end
end

acceleration_label = @sprintf("%02dg", round(Int, 10 * acceleration_rms_g))
output_path = get(ENV, "CANTILEVER_SWEEP_OUTPUT",
                  joinpath(@__DIR__, "reference",
                           "base_excited_cantilever_sweep_" *
                           acceleration_label * ".csv"))
mkpath(dirname(output_path))
open(output_path, "w") do io
    println(io,
            "acceleration_rms_g,frequency_hz,frequency_over_discrete_f1," *
            "omega,transverse_peak,longitudinal_minimum,rotation_peak," *
            "periodicity_error,cycles,accepted_steps,rejected_steps")
    for row in rows
        @printf(io,
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d\n",
                row.acceleration_rms_g, row.frequency_hz,
                row.frequency_over_discrete_f1, row.omega,
                row.transverse_peak, row.longitudinal_minimum,
                row.rotation_peak, row.periodicity_error, row.cycles,
                row.accepted_steps, row.rejected_steps)
    end
end

metadata_path = replace(output_path, ".csv" => "_metadata.txt")
open(metadata_path, "w") do io
    println(io, "generated_at = ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "polydeg = ", polydeg)
    println(io, "refinement_level = ", refinement_level)
    println(io, "cells = ", 2^refinement_level)
    println(io, "constraint_multiplier = ", constraint_multiplier)
    println(io, "linear_eigenvalue = ", mode)
    println(io, "linear_frequency_hz = ", linear_frequency_hz)
    println(io, "linear_damping_ratio = ", linear_damping_ratio)
    println(io, "relative_tolerance = ",
            get(ENV, "CANTILEVER_SWEEP_RELTOL", "1e-4"))
    println(io, "absolute_tolerance = ",
            get(ENV, "CANTILEVER_SWEEP_ABSTOL", "1e-6"))
    println(io, "steps_per_period_limit = ",
            get(ENV, "CANTILEVER_STEPS_PER_PERIOD", "32"))
    println(io, "source = Farokhi, Xia, and Erturk (2022), DOI 10.1007/s11071-021-07023-9")
end

println("Wrote ", output_path)
println("Wrote ", metadata_path)
