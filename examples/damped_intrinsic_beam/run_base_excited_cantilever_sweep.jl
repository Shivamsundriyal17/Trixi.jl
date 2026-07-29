using Dates
using LinearAlgebra: eigvals
using Printf
using Serialization: deserialize, serialize

# Load the benchmark once. The short setup solve keeps the elixir directly
# executable while avoiding an additional long fixed-frequency calculation
# before the continuation sweep.
ENV["CANTILEVER_SETTLING_CYCLES"] = "1"
ENV["CANTILEVER_MEASUREMENT_CYCLES"] = "2"
ENV["CANTILEVER_RAMP_CYCLES"] = "0"
ENV["CANTILEVER_SETUP_ONLY"] = "true"
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

function final_cycle_ledger(solution)
    first_index = length(solution.u) - samples_per_cycle
    return cantilever_cycle_ledger(solution.u[first_index:end],
                                   solution.t[first_index:end])
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
reference_frequency_hz = FAROKHI_REFERENCE_OMEGA1 /
                         (2.0 * pi * FAROKHI_T)
frequencies = normalized_frequencies .* linear_frequency_hz
@printf("Unloaded Euler-Bernoulli reference frequency = %.8f Hz\n",
        reference_frequency_hz)
@printf("Discrete gravity-loaded mode: lambda = %.8e %+.8ei, f1 = %.8f Hz, zeta = %.6e\n",
        real(mode), imag(mode), linear_frequency_hz, linear_damping_ratio)

acceleration_label = @sprintf("%02dg", round(Int, 10 * acceleration_rms_g))
output_path = get(ENV, "CANTILEVER_SWEEP_OUTPUT",
                  joinpath(@__DIR__, "reference",
                           "base_excited_cantilever_sweep_" *
                           acceleration_label * ".csv"))
checkpoint_path = get(ENV, "CANTILEVER_SWEEP_CHECKPOINT",
                      replace(output_path, ".csv" => "_checkpoint.jls"))
mkpath(dirname(output_path))
mkpath(dirname(checkpoint_path))

function write_rows(path, rows)
    temporary_path = path * ".tmp"
    open(temporary_path, "w") do io
        println(io,
                "acceleration_rms_g,frequency_hz,normalized_frequency," *
                "frequency_over_unloaded_omega1,omega,transverse_peak," *
                "longitudinal_minimum,rotation_peak," *
                "periodicity_error,cycles,accepted_steps,rejected_steps," *
                "total_energy_change,material_dissipation,jump_dissipation," *
                "left_boundary_dissipation,right_boundary_dissipation," *
                "physical_root_work,sat_data_work,ledger_residual," *
                "relative_ledger_residual")
        for row in rows
            @printf(io,
                    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    row.acceleration_rms_g, row.frequency_hz,
                    row.normalized_frequency,
                    row.frequency_over_unloaded_omega1, row.omega,
                    row.transverse_peak, row.longitudinal_minimum,
                    row.rotation_peak, row.periodicity_error, row.cycles,
                    row.accepted_steps, row.rejected_steps,
                    row.total_energy_change, row.material_dissipation,
                    row.jump_dissipation, row.left_boundary_dissipation,
                    row.right_boundary_dissipation, row.physical_root_work,
                    row.sat_data_work, row.ledger_residual,
                    row.relative_ledger_residual)
        end
    end
    mv(temporary_path, path; force = true)
    return nothing
end

function write_checkpoint(path, rows, continuation_state, next_index;
                          partial_cycles = 0)
    temporary_path = path * ".tmp"
    checkpoint = (;
                  format_version = 4,
                  acceleration_rms_g,
                  polydeg,
                  refinement_level,
                  constraint_multiplier,
                  normalized_frequencies,
                  reference_frequency_hz,
                  linear_frequency_hz,
                  rows,
                  continuation_state,
                  next_index,
                  partial_cycles)
    serialize(temporary_path, checkpoint)
    mv(temporary_path, path; force = true)

    archive_states = lowercase(get(ENV,
                                   "CANTILEVER_SWEEP_ARCHIVE_STATES",
                                   "true")) in ("1", "true", "yes")
    if archive_states && iszero(partial_cycles)
        stem, extension = splitext(path)
        point_path = @sprintf("%s_point_%03d%s", stem, next_index - 1,
                              extension)
        point_temporary_path = point_path * ".tmp"
        serialize(point_temporary_path, checkpoint)
        mv(point_temporary_path, point_path; force = true)
    end
    return nothing
end

resume_requested = lowercase(get(ENV, "CANTILEVER_SWEEP_RESUME",
                                 "false")) in ("1", "true", "yes")
rows = NamedTuple[]
continuation_state = copy(split_problem.u0)
first_index = 1
partial_cycles = 0
if resume_requested && isfile(checkpoint_path)
    checkpoint = deserialize(checkpoint_path)
    checkpoint.format_version in (1, 3, 4) ||
        error("unsupported sweep checkpoint format")
    checkpoint.acceleration_rms_g == acceleration_rms_g ||
        error("checkpoint acceleration does not match this run")
    checkpoint.polydeg == polydeg ||
        error("checkpoint polynomial degree does not match this run")
    checkpoint.refinement_level == refinement_level ||
        error("checkpoint refinement level does not match this run")
    checkpoint.constraint_multiplier == constraint_multiplier ||
        error("checkpoint constraint multiplier does not match this run")
    checkpoint.normalized_frequencies == normalized_frequencies ||
        error("checkpoint frequency path does not match this run")
    checkpoint.linear_frequency_hz == linear_frequency_hz ||
        error("checkpoint discrete frequency does not match this run")
    if checkpoint.format_version in (3, 4)
        checkpoint.reference_frequency_hz == reference_frequency_hz ||
            error("checkpoint reference frequency does not match this run")
        rows = checkpoint.rows
    else
        # Version 1 already used the correct gravity-loaded normalization,
        # but did not retain the unloaded Euler-Bernoulli diagnostic column.
        rows = NamedTuple[merge(row,
                                (normalized_frequency =
                                     row.frequency_over_discrete_f1,
                                 frequency_over_unloaded_omega1 =
                                     row.frequency_hz /
                                     reference_frequency_hz))
                          for row in checkpoint.rows]
    end
    continuation_state = checkpoint.continuation_state
    first_index = checkpoint.next_index
    partial_cycles = checkpoint.format_version == 4 ?
                     checkpoint.partial_cycles : 0
    @printf("Resuming %s at point %d of %d after %d partial cycles\n",
            checkpoint_path, first_index, length(frequencies),
            partial_cycles)
end

let continuation_state = continuation_state,
    resumed_partial_cycles = partial_cycles
    for index in first_index:length(frequencies)
        frequency = frequencies[index]
        is_partial_resume = index == first_index &&
                            resumed_partial_cycles > 0
        cycles = is_partial_resume ? additional_cycles :
                 index == 1 ? first_cycles : continuation_cycles
        ramp = is_partial_resume ? 0.0 :
               index == 1 ?
               parse(Float64,
                     get(ENV, "CANTILEVER_SWEEP_RAMP_CYCLES", "8")) :
               0.0
        local_solution, metrics = solve_cycles(continuation_state, cycles,
                                               frequency;
                                               ramp_cycles = ramp)
        continuation_state = copy(local_solution.u[end])
        local total_cycles = resumed_partial_cycles + cycles
        resumed_partial_cycles = 0
        @printf("  point %d/%d settling: %d cycles, periodicity %.3e\n",
                index, length(frequencies), total_cycles,
                metrics.periodicity_error)
        flush(stdout)

        while metrics.periodicity_error > periodicity_tolerance &&
              total_cycles < maximum_cycles
            # Save a same-frequency restart before every additional block.
            # If the next block is interrupted, only that block is repeated.
            write_checkpoint(checkpoint_path, rows, continuation_state,
                             index; partial_cycles = total_cycles)
            extra_solution, metrics = solve_cycles(continuation_state,
                                                    additional_cycles,
                                                    frequency)
            continuation_state = copy(extra_solution.u[end])
            total_cycles += additional_cycles
            local_solution = extra_solution
            @printf("  point %d/%d settling: %d cycles, periodicity %.3e\n",
                    index, length(frequencies), total_cycles,
                    metrics.periodicity_error)
            flush(stdout)
        end

        ledger = final_cycle_ledger(local_solution)
        row = (;
               acceleration_rms_g,
               frequency_hz = frequency,
               normalized_frequency = frequency / linear_frequency_hz,
               frequency_over_unloaded_omega1 =
                   frequency / reference_frequency_hz,
               omega = root_velocity.omega,
               transverse_peak = metrics.transverse_peak,
               longitudinal_minimum = metrics.longitudinal_minimum,
               rotation_peak = metrics.rotation_peak,
               periodicity_error = metrics.periodicity_error,
               cycles = total_cycles,
               accepted_steps = local_solution.destats.naccept,
               rejected_steps = local_solution.destats.nreject,
               ledger...)
        push!(rows, row)
        write_rows(output_path, rows)
        write_checkpoint(checkpoint_path, rows, continuation_state,
                         index + 1; partial_cycles = 0)
        @printf("f=%7.4f Hz  Omega/omega1=%8.5f  f/f1_unloaded=%8.5f  |w|=%9.6f  u_min=%9.6f  |psi|=%9.6f  periodic=%8.2e  ledger=%8.2e  cycles=%d\n",
                row.frequency_hz, row.normalized_frequency,
                row.frequency_over_unloaded_omega1,
                row.transverse_peak, row.longitudinal_minimum,
                row.rotation_peak, row.periodicity_error,
                row.relative_ledger_residual, row.cycles)
    end
end

write_rows(output_path, rows)

metadata_path = replace(output_path, ".csv" => "_metadata.txt")
open(metadata_path, "w") do io
    println(io, "generated_at = ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "polydeg = ", polydeg)
    println(io, "refinement_level = ", refinement_level)
    println(io, "cells = ", 2^refinement_level)
    println(io, "constraint_multiplier = ", constraint_multiplier)
    println(io, "linear_eigenvalue = ", mode)
    println(io, "reference_frequency_hz = ", reference_frequency_hz)
    println(io, "linear_frequency_hz = ", linear_frequency_hz)
    println(io, "linear_damping_ratio = ", linear_damping_ratio)
    println(io, "normalized_frequency_definition = frequency / gravity-loaded discrete first mode")
    println(io, "normalized_frequencies = ", join(normalized_frequencies, ","))
    println(io, "first_cycles = ", first_cycles)
    println(io, "continuation_cycles = ", continuation_cycles)
    println(io, "additional_cycles = ", additional_cycles)
    println(io, "maximum_cycles = ", maximum_cycles)
    println(io, "periodicity_tolerance = ", periodicity_tolerance)
    println(io, "archive_states = ",
            get(ENV, "CANTILEVER_SWEEP_ARCHIVE_STATES", "true"))
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
