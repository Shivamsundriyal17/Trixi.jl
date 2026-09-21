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
const SUMMARY_ROW_FORMAT = Format("%s,%.1f,%s,%.1f,%d,%d,%.4f," *
                                  "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e," *
                                  "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e," *
                                  "%.16e,%.16e,%.16e,%.16e,%.16e\n")

output_directory = isempty(ARGS) ?
                   joinpath(EXAMPLE_DIRECTORY, "results",
                            "rotating_beam_campaign") :
                   abspath(ARGS[1])
summary_file = length(ARGS) < 2 ?
               joinpath(output_directory, "summary.csv") :
               abspath(ARGS[2])
mkpath(output_directory)

quick_mode = lowercase(get(ENV, "ROTATING_CAMPAIGN_QUICK", "false")) in ("1", "true", "yes")
reuse_results = lowercase(get(ENV, "ROTATING_REUSE_RESULTS", "false")) in ("1", "true",
                                                                           "yes")
transient_final_time = quick_mode ? 1.0e-4 : 40.0
steady_final_time = quick_mode ? 1.0e-4 : 10.0
save_count = quick_mode ? 2 : 401

cases = [
    (name = "baseline", damping = 1.0, steady = false,
     final_time = transient_final_time,
     record_ledger = true),
    (name = "double_damping", damping = 2.0, steady = false,
     final_time = transient_final_time, record_ledger = false),
    (name = "undamped", damping = 0.0, steady = false,
     final_time = transient_final_time,
     record_ledger = false),
    (name = "steady_initial", damping = 1.0, steady = true,
     final_time = steady_final_time, record_ledger = true)
]

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
for case in cases
    output_file = joinpath(output_directory, case.name * ".jls")
    configuration = (; case_name = case.name, damping_multiplier = case.damping,
                     steady_initial_condition = case.steady,
                     polydeg = parse(Int, get(ENV, "ROTATING_POLYDEG", "3")),
                     ncells = 2^parse(Int, get(ENV, "ROTATING_REFINEMENT_LEVEL", "3")),
                     cfl = parse(Float64, get(ENV, "ROTATING_CFL", "0.01")),
                     final_time = case.final_time, save_count,
                     record_online_ledger = case.record_ledger)
    if reuse_results && isfile(output_file)
        println("reusing ", output_file)
    else
        @printf("rotating case=%s, damping=%g, steady_initial=%s, T=%g\n",
                case.name, case.damping, case.steady, case.final_time)
        withenv("ROTATING_CASE" => case.name,
                "ROTATING_DAMPING_MULTIPLIER" => string(case.damping),
                "ROTATING_STEADY_INITIAL" => string(case.steady),
                "ROTATING_T_END" => string(case.final_time),
                "ROTATING_SAVE_COUNT" => string(save_count),
                "ROTATING_RECORD_LEDGER" => string(case.record_ledger),
                "ROTATING_OUTPUT_FILE" => output_file) do
            case_module = Module(Symbol("RotatingBeam_", case.name))
            Base.include(case_module, EXPORTER)
        end
    end
    push!(results,
          validate_rotating_result(deserialize(output_file),
                                   configuration, provenance))
    GC.gc()
end

commit, worktree_state = repository_metadata()
mkpath(dirname(summary_file))
open(summary_file, "w") do io
    println(io, "# Clean Trixi rotating-beam campaign")
    println(io, "# generated_utc=", now(UTC))
    println(io, "# julia_version=", VERSION)
    println(io, "# trixi_commit=", commit)
    println(io, "# trixi_worktree=", worktree_state)
    println(io, "# julia_threads=", Threads.nthreads())
    println(io, "# hyperbolic_flux=characteristic_upwind")
    println(io, "# auxiliary_flux=alternating_ldg")
    println(io,
            "case,damping_multiplier,steady_initial,T,polydeg,ncells,cfl," *
            "energy_final,energy_steady,energy_relative_error," *
            "f1_relative_l2,f1_relative_linf,v2_relative_l2," *
            "v2_relative_linf,late_f1_relative_linf," *
            "late_v2_relative_linf,beam_length_relative_error," *
            "maximum_ledger_relative_residual," *
            "cumulative_material_dissipation," *
            "cumulative_jump_dissipation," *
            "cumulative_left_boundary_dissipation," *
            "cumulative_right_boundary_dissipation," *
            "cumulative_physical_root_work," *
            "cumulative_sat_data_work," *
            "cumulative_ledger_residual," *
            "cumulative_ledger_relative_residual")
    for data in results
        diagnostics = data.diagnostics
        metrics = diagnostics.final_metrics
        cumulative_recorded = isfinite(diagnostics.cumulative_ledger_residual)
        cumulative_terms = cumulative_recorded ?
                           (diagnostics.cumulative_material_dissipation,
                            diagnostics.cumulative_jump_dissipation,
                            diagnostics.cumulative_left_boundary_dissipation,
                            diagnostics.cumulative_right_boundary_dissipation,
                            diagnostics.cumulative_physical_root_work,
                            diagnostics.cumulative_sat_data_work) :
                           ntuple(_ -> NaN, 6)
        cumulative_relative_residual = if cumulative_recorded
            cumulative_scale = max(abs(diagnostics.energy_final -
                                       first(data.energy_history)),
                                   sum(abs, cumulative_terms), 1.0)
            abs(diagnostics.cumulative_ledger_residual) / cumulative_scale
        else
            NaN
        end
        format(io, SUMMARY_ROW_FORMAT,
               data.case_name, data.damping_multiplier,
               data.steady_initial_condition, last(data.times),
               data.polydeg, data.ncells, data.cfl,
               diagnostics.energy_final, diagnostics.energy_steady,
               diagnostics.energy_relative_error,
               metrics.f1_relative_l2, metrics.f1_relative_linf,
               metrics.v2_relative_l2, metrics.v2_relative_linf,
               diagnostics.late_f1_relative_linf,
               diagnostics.late_v2_relative_linf,
               diagnostics.beam_length_relative_error,
               diagnostics.maximum_ledger_relative_residual,
               cumulative_terms...,
               diagnostics.cumulative_ledger_residual,
               cumulative_relative_residual)
    end
end

println("wrote rotating-beam campaign to ", output_directory)
println("wrote ", summary_file)
