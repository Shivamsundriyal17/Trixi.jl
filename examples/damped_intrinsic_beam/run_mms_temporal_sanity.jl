using Dates: UTC, now
using Printf: @printf
using Trixi

const ELIXIR = joinpath(@__DIR__, "elixir_mms_physical.jl")
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_list(name, default, conversion)
    values = split(get(ENV, name, default), ',')
    return conversion.(strip.(values))
end

polydeg = parse(Int, get(ENV, "MMS_TEMPORAL_POLYDEG", "3"))
refinement_level = parse(Int,
                         get(ENV, "MMS_TEMPORAL_REFINEMENT_LEVEL", "7"))
sigma = parse(Float64, get(ENV, "MMS_TEMPORAL_SIGMA", "1.0"))
auxiliary_flux = get(ENV, "MMS_TEMPORAL_AUXILIARY_FLUX", "alternating")
tolerances = parse_list("MMS_TEMPORAL_TOLS", "1e-8,1e-10,1e-12,1e-13",
                        value -> parse(Float64, value))
output_file = isempty(ARGS) ?
              joinpath(@__DIR__, "results", "mms_temporal_sanity.csv") :
              abspath(ARGS[1])

module MMSTemporalRun
using Trixi
end

function parabolic_solver(name)
    name == "alternating" && return ViscousFormulationLocalDG()
    name == "br1" && return ViscousFormulationBassiRebay1()
    throw(ArgumentError("unknown auxiliary flux '$name'; use alternating or br1"))
end

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

function run_case(tolerance)
    Trixi.trixi_include(MMSTemporalRun, ELIXIR;
                        polydeg,
                        initial_refinement_level = refinement_level,
                        sigma,
                        solver_parabolic = parabolic_solver(auxiliary_flux),
                        time_int_tol = tolerance)
    errors = Base.invokelatest(MMSTemporalRun.analysis_callback,
                               MMSTemporalRun.sol)
    l2 = collect(errors.l2)
    linf = collect(errors.linf)
    return (u1_l2 = sum(l2[1:6]) / 6,
            u2_l2 = sum(l2[7:12]) / 6,
            maximum_component_l2 = maximum(l2),
            maximum_component_linf = maximum(linf))
end

mkpath(dirname(output_file))
commit, worktree_state = repository_metadata()
ncells = 2^refinement_level

open(output_file, "w") do io
    println(io, "# Physical MMS temporal-error sanity check")
    println(io, "# generated_utc=", now(UTC))
    println(io, "# julia_version=", VERSION)
    println(io, "# trixi_commit=", commit)
    println(io, "# trixi_worktree=", worktree_state)
    println(io, "# julia_threads=", Threads.nthreads())
    println(io, "# polydeg=", polydeg)
    println(io, "# refinement_level=", refinement_level)
    println(io, "# ncells=", ncells)
    println(io, "# sigma=", sigma)
    println(io, "# auxiliary_flux=", auxiliary_flux)
    println(io,
            "tolerance,u1_mean_l2,u2_mean_l2,maximum_component_l2," *
            "maximum_component_linf")

    for tolerance in tolerances
        @printf("MMS temporal sanity: k=%d, N=%d, tolerance=%.1e ... ",
                polydeg, ncells, tolerance)
        flush(stdout)
        errors = run_case(tolerance)
        @printf("u1=%.8e, u2=%.8e\n", errors.u1_l2, errors.u2_l2)
        @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e\n",
                tolerance, errors.u1_l2, errors.u2_l2,
                errors.maximum_component_l2,
                errors.maximum_component_linf)
        flush(io)
        GC.gc()
    end
end

println("wrote ", output_file)
