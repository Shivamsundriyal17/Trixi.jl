module BeamRunHelpers

using Trixi
using SciMLBase: successful_retcode
using LinearAlgebra: opnorm
using SHA: sha256

export require_complete_solution, checkpoint_gravity, require_matching_configuration,
       beam_explicit_timestep, beam_source_identity, beam_run_provenance,
       validate_rotating_result, beam_node, cantilever_checkpoint_environment

"""Reject failed, partial, or nonfinite solutions before reporting numerical results."""
function require_complete_solution(sol, final_time)
    successful_retcode(sol.retcode) ||
        error("beam solve failed with return code $(sol.retcode)")
    isempty(sol.t) && error("beam solve returned no saved times")
    sol.t[end] == final_time ||
        error("beam solve ended at $(sol.t[end]); expected $final_time")
    all(state -> all(isfinite, state), sol.u) ||
        error("beam solve contains nonfinite saved states")
    return sol
end

"""Older checkpoints predate configurable gravity and always used multiplier one."""
function checkpoint_gravity(checkpoint)
    value = hasproperty(checkpoint, :gravity_multiplier) ?
            checkpoint.gravity_multiplier : 1.0
    isfinite(value) && value >= 0 || error("invalid checkpoint gravity multiplier")
    return value
end

"""Restore physical and (for format 5) numerical settings for checkpoint replay."""
function cantilever_checkpoint_environment(checkpoint, row)
    environment = Dict("CANTILEVER_ACCELERATION_RMS_G" => string(checkpoint.acceleration_rms_g),
                       "CANTILEVER_FREQUENCY_HZ" => string(row.frequency_hz),
                       "CANTILEVER_POLYDEG" => string(checkpoint.polydeg),
                       "CANTILEVER_REFINEMENT_LEVEL" => string(checkpoint.refinement_level),
                       "CANTILEVER_CONSTRAINT_MULTIPLIER" => string(checkpoint.constraint_multiplier),
                       "CANTILEVER_GRAVITY_MULTIPLIER" => string(checkpoint_gravity(checkpoint)))
    if hasproperty(checkpoint, :run_configuration)
        config = checkpoint.run_configuration
        environment["CANTILEVER_RELTOL"] = string(config.reltol)
        environment["CANTILEVER_ABSTOL"] = string(config.abstol)
        environment["CANTILEVER_DTMAX"] = string(2pi / row.omega / config.steps_per_period)
        environment["CANTILEVER_SPARSE_JACOBIAN"] = string(config.use_sparse_jacobian)
    end
    return environment
end

function require_matching_configuration(saved, requested)
    for key in keys(requested)
        hasproperty(saved, key) || error("saved run lacks configuration field $key")
        isequal(getproperty(saved, key), getproperty(requested, key)) ||
            error("saved run configuration differs for $key")
    end
    return nothing
end

"""
    beam_explicit_timestep(ode, wave_cfl; diffusion_cfl=0.1)

Limit the wave CFL step by a conservative principal-diffusion estimate
`diffusion_cfl * h_min^2 / (norm(M^-1 H, Inf) * (p+1)^4)`.
This is a fixed-mesh explicit-step safeguard, not a nonlinear stability theorem.
The default is checked against assembled beam operators in the beam regression suite.
"""
function beam_explicit_timestep(ode, wave_cfl; diffusion_cfl = 0.1)
    isfinite(wave_cfl) && wave_cfl > 0 ||
        throw(ArgumentError("wave CFL must be positive and finite"))
    isfinite(diffusion_cfl) && 0 < diffusion_cfl <= 0.1 ||
        throw(ArgumentError("diffusion CFL must be in (0, 0.1]"))
    semi = ode.p
    equations = semi.equations
    diffusivity = opnorm(equations.mass_inverse * equations.damping_operator, Inf)
    wave_step = StepsizeCallback(cfl = wave_cfl)(ode)
    iszero(diffusivity) && return wave_step
    h_min = 2 / maximum(abs, semi.cache.elements.inverse_jacobian)
    diffusion_step = diffusion_cfl * h_min^2 /
                     (diffusivity * Trixi.nnodes(semi.solver)^4)
    return min(wave_step, diffusion_step)
end

@inline function beam_node(array, i, element)
    return SVector{12}(ntuple(variable -> array[variable, i, element], Val(12)))
end

"""Fingerprint executable beam sources and pinned environments, excluding results."""
function beam_source_identity(root = normpath(joinpath(@__DIR__, "..", "..")))
    paths = String[]
    for relative in ("src", "examples/damped_intrinsic_beam")
        for (directory, subdirectories, files) in walkdir(joinpath(root, relative))
            filter!(name -> !(name in ("results", "reference", "visualization",
                                       "__pycache__")), subdirectories)
            for name in files
                (endswith(name, ".jl") || endswith(name, ".toml")) &&
                    push!(paths, relpath(joinpath(directory, name), root))
            end
        end
    end
    io = IOBuffer()
    for path in sort!(paths)
        write(io, path, '\0', read(joinpath(root, path)), '\0')
    end
    return bytes2hex(sha256(take!(io)))
end

function beam_run_provenance()
    root = normpath(joinpath(@__DIR__, "..", ".."))
    commit = try
        readchomp(`git -C $root rev-parse HEAD`)
    catch
        "unavailable"
    end
    return (; source_identity = beam_source_identity(), julia_version = string(VERSION),
            threads = Threads.nthreads(), generating_commit = commit)
end

function validate_rotating_result(data, configuration, provenance)
    hasproperty(data, :run_configuration) ||
        error("legacy rotating results lack validation metadata; regenerate them")
    require_matching_configuration(data.run_configuration, configuration)
    require_matching_configuration(data.provenance,
                                   (; source_identity = provenance.source_identity,
                                    julia_version = provenance.julia_version,
                                    threads = provenance.threads))
    data.completed || error("cached rotating solve did not complete")
    !isempty(data.times) && data.times[end] == configuration.final_time ||
        error("cached rotating solve has the wrong final time")
    length(data.times) == length(data.states) == configuration.save_count ||
        error("cached rotating solve has the wrong number of saved states")
    all(state -> all(isfinite, state), data.states) ||
        error("cached rotating states contain nonfinite values")
    return data
end

end # module
