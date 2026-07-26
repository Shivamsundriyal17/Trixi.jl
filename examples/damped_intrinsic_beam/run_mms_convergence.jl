using Dates: now, UTC
using Printf: @printf, @sprintf
using Trixi

const ELIXIR = joinpath(@__DIR__, "elixir_mms_physical.jl")
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_list(name, default, conversion)
    values = split(get(ENV, name, default), ',')
    return conversion.(strip.(values))
end

polydegs = parse_list("MMS_POLYDEGS", "1,2,3", value -> parse(Int, value))
refinement_levels = parse_list("MMS_REFINEMENT_LEVELS", "2,3,4,5,6,7",
                               value -> parse(Int, value))
penalties = parse_list("MMS_SIGMAS", "1", value -> parse(Float64, value))
auxiliary_fluxes = parse_list("MMS_AUXILIARY_FLUXES", "alternating,br1",
                              identity)
time_int_tol = parse(Float64, get(ENV, "MMS_TIME_TOL", "1e-12"))
output_file = isempty(ARGS) ?
              joinpath(@__DIR__, "results", "mms_convergence.csv") :
              abspath(ARGS[1])

module MMSRun
using Trixi
end

function parabolic_solver(name)
    name == "alternating" && return ViscousFormulationLocalDG()
    name == "br1" && return ViscousFormulationBassiRebay1()
    throw(ArgumentError("unknown auxiliary flux '$name'; use alternating or br1"))
end

function hyperbolic_flux_name(sigma)
    sigma == 1 && return "upwind"
    sigma == 0 && return "central"
    return @sprintf("sigma_%g", sigma)
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

function run_case(polydeg, refinement_level, sigma, auxiliary_flux, tolerance)
    Trixi.trixi_include(MMSRun, ELIXIR;
                        polydeg,
                        initial_refinement_level = refinement_level,
                        sigma,
                        solver_parabolic = parabolic_solver(auxiliary_flux),
                        time_int_tol = tolerance)
    errors = Base.invokelatest(MMSRun.analysis_callback, MMSRun.sol)
    error_u1 = sum(errors.l2[1:6]) / 6
    error_u2 = sum(errors.l2[7:12]) / 6
    return error_u1, error_u2
end

mkpath(dirname(output_file))
commit, worktree_state = repository_metadata()
previous = Dict{Tuple{Int, Float64, String, String}, Tuple{Int, Float64}}()

open(output_file, "w") do io
    println(io, "# Physical, constitutively consistent Kelvin-Voigt MMS")
    println(io, "# generated_utc=", now(UTC))
    println(io, "# julia_version=", VERSION)
    println(io, "# trixi_commit=", commit)
    println(io, "# trixi_worktree=", worktree_state)
    println(io, "# julia_threads=", Threads.nthreads())
    println(io, "# command=", join(Base.ARGS, ' '))
    println(io, "# polydegs=", join(polydegs, ','))
    println(io, "# refinement_levels=", join(refinement_levels, ','))
    println(io, "# sigmas=", join(penalties, ','))
    println(io, "# auxiliary_fluxes=", join(auxiliary_fluxes, ','))
    println(io, "# rock4_abstol_reltol=", time_int_tol)
    println(io, "# error=arithmetic mean of six componentwise L2 errors at T=1")
    println(io,
            "polydeg,refinement_level,ncells,h,hyperbolic_flux,auxiliary_flux," *
            "block,l2_error,eoc")

    for polydeg in polydegs, sigma in penalties,
        auxiliary_flux in auxiliary_fluxes,
        refinement_level in refinement_levels
        ncells = 2^refinement_level
        h = 1 / ncells
        flux_name = hyperbolic_flux_name(sigma)
        @printf("k=%d, N=%d, hyperbolic=%s, auxiliary=%s ... ",
                polydeg, ncells, flux_name, auxiliary_flux)
        flush(stdout)

        error_u1, error_u2 = run_case(polydeg, refinement_level, sigma,
                                      auxiliary_flux, time_int_tol)
        @printf("u1=%.8e, u2=%.8e\n", error_u1, error_u2)

        for (block, error) in (("u1", error_u1), ("u2", error_u2))
            key = (polydeg, sigma, auxiliary_flux, block)
            eoc = if haskey(previous, key)
                previous_cells, previous_error = previous[key]
                log(previous_error / error) / log(ncells / previous_cells)
            else
                NaN
            end
            @printf(io, "%d,%d,%d,%.16e,%s,%s,%s,%.16e,%.8f\n",
                    polydeg, refinement_level, ncells, h, flux_name,
                    auxiliary_flux, block, error, eoc)
            previous[key] = (ncells, error)
        end
        flush(io)
        GC.gc()
    end
end

println("wrote ", output_file)
