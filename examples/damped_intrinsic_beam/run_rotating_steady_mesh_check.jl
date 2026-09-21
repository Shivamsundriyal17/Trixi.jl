if !isdefined(@__MODULE__, :BeamRunHelpers)
    Base.include(@__MODULE__, joinpath(@__DIR__, "beam_run_helpers.jl"))
end
using .BeamRunHelpers

using Dates: UTC, now
using Printf: Format, @printf, format
using Serialization: deserialize

const EXAMPLE_DIRECTORY = @__DIR__
const EXPORTER = joinpath(EXAMPLE_DIRECTORY, "save_rotating_beam_data.jl")
const REPOSITORY_ROOT = normpath(joinpath(EXAMPLE_DIRECTORY, "..", ".."))
const SUMMARY_ROW_FORMAT = Format("%d,%d,%.4f,%.16e,%.16e,%.16e,%.16e,%.16e," *
                                  "%.16e,%.16e,%.8f,%.8f\n")

function parse_list(name, default, conversion)
    values = split(get(ENV, name, default), ',')
    return conversion.(strip.(values))
end

refinement_levels = parse_list("ROTATING_STEADY_MESH_LEVELS", "2,3,4",
                               value -> parse(Int, value))
final_time = parse(Float64,
                   get(ENV, "ROTATING_STEADY_MESH_T_END", "1.0"))
reuse_results = lowercase(get(ENV, "ROTATING_STEADY_REUSE_RESULTS",
"false")) in ("1", "true", "yes")
output_directory = isempty(ARGS) ?
                   joinpath(EXAMPLE_DIRECTORY, "results",
                            "rotating_steady_mesh_check") :
                   abspath(ARGS[1])
summary_file = length(ARGS) < 2 ?
               joinpath(output_directory, "summary.csv") :
               abspath(ARGS[2])
mkpath(output_directory)

function repository_metadata()
    commit = try
        readchomp(`git -C $REPOSITORY_ROOT rev-parse HEAD`)
    catch
        "unavailable"
    end
    status = try
        readchomp(`git -C $REPOSITORY_ROOT status --porcelain`)
    catch
        "unavailable"
    end
    state = status == "unavailable" ? status :
            (isempty(status) ? "clean" : "dirty")
    return commit, state
end

provenance = beam_run_provenance()
results = []
for refinement_level in refinement_levels
    ncells = 2^refinement_level
    output_file = joinpath(output_directory, "steady_N$(ncells).jls")
    configuration = (; case_name = "steady_N$(ncells)", damping_multiplier = 1.0,
                     steady_initial_condition = true,
                     polydeg = parse(Int, get(ENV, "ROTATING_POLYDEG", "3")), ncells,
                     cfl = parse(Float64, get(ENV, "ROTATING_CFL", "0.01")),
                     final_time, save_count = 2, record_online_ledger = false)
    if reuse_results && isfile(output_file)
        @printf("reusing rotating steady mesh check N=%d\n", ncells)
    else
        @printf("rotating steady mesh check: N=%d, T=%g\n",
                ncells, final_time)
        withenv("ROTATING_CASE" => "steady_N$(ncells)",
                "ROTATING_DAMPING_MULTIPLIER" => "1.0",
                "ROTATING_STEADY_INITIAL" => "true",
                "ROTATING_REFINEMENT_LEVEL" => string(refinement_level),
                "ROTATING_RECORD_LEDGER" => "false",
                "ROTATING_T_END" => string(final_time),
                "ROTATING_SAVE_COUNT" => "2",
                "ROTATING_OUTPUT_FILE" => output_file) do
            case_module = Module(Symbol("RotatingSteady_N", ncells))
            Base.include(case_module, EXPORTER)
        end
    end
    push!(results,
          validate_rotating_result(deserialize(output_file),
                                   configuration, provenance))
    GC.gc()
end

commit, worktree_state = repository_metadata()
open(summary_file, "w") do io
    println(io, "# Analytic rotating-branch steady-initialized mesh check")
    println(io, "# generated_utc=", now(UTC))
    println(io, "# julia_version=", VERSION)
    println(io, "# trixi_commit=", commit)
    println(io, "# trixi_worktree=", worktree_state)
    println(io, "# julia_threads=", Threads.nthreads())
    println(io, "# final_time=", final_time)
    println(io, "# hyperbolic_flux=characteristic_upwind")
    println(io, "# auxiliary_flux=alternating_ldg")
    println(io,
            "ncells,polydeg,cfl,f1_relative_l2,f1_relative_linf," *
            "v2_relative_l2,v2_relative_linf,energy_relative_error," *
            "beam_length_relative_error,maximum_ledger_relative_residual," *
            "f1_linf_eoc,v2_linf_eoc")

    for (index, data) in enumerate(results)
        metrics = data.diagnostics.final_metrics
        previous = index == 1 ? nothing : results[index - 1]
        previous_metrics = isnothing(previous) ?
                           nothing : previous.diagnostics.final_metrics
        f1_eoc = isnothing(previous) ? NaN :
                 log(previous_metrics.f1_relative_linf /
                     metrics.f1_relative_linf) /
                 log(data.ncells / previous.ncells)
        v2_eoc = isnothing(previous) ? NaN :
                 log(previous_metrics.v2_relative_linf /
                     metrics.v2_relative_linf) /
                 log(data.ncells / previous.ncells)
        format(io, SUMMARY_ROW_FORMAT,
               data.ncells, data.polydeg, data.cfl,
               metrics.f1_relative_l2, metrics.f1_relative_linf,
               metrics.v2_relative_l2, metrics.v2_relative_linf,
               data.diagnostics.energy_relative_error,
               data.diagnostics.beam_length_relative_error,
               data.diagnostics.maximum_ledger_relative_residual,
               f1_eoc, v2_eoc)
    end
end

println("wrote rotating steady mesh check to ", summary_file)
