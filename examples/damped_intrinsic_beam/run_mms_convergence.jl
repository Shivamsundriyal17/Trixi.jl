using Dates: now, UTC
using Printf: Format, @printf, @sprintf, format
using Trixi

const ELIXIR = joinpath(@__DIR__, "elixir_mms_physical.jl")
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const COMPONENT_ROW_FORMAT = Format("%d,%d,%d,%.16e,%s,%s,%d,%s,%s," *
                                    "%.16e,%.8f,%.16e,%.8f\n")

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
component_output_file = length(ARGS) < 2 ?
                        first(splitext(output_file)) *
                        "_componentwise.csv" :
                        abspath(ARGS[2])
component_names = ("v1", "v2", "v3", "omega1", "omega2", "omega3",
                   "f1", "f2", "f3", "m1", "m2", "m3")

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
    return (l2 = collect(errors.l2), linf = collect(errors.linf))
end

mkpath(dirname(output_file))
mkpath(dirname(component_output_file))
commit, worktree_state = repository_metadata()
previous = Dict{Tuple{Int, Float64, String, String}, Tuple{Int, Float64}}()
previous_component_l2 = Dict{Tuple{Int, Float64, String, Int},
                             Tuple{Int, Float64}}()
previous_component_linf = Dict{Tuple{Int, Float64, String, Int},
                               Tuple{Int, Float64}}()

component_io = open(component_output_file, "w")
try
    println(component_io,
            "# Componentwise physical Kelvin-Voigt MMS errors")
    println(component_io, "# generated_utc=", now(UTC))
    println(component_io, "# julia_version=", VERSION)
    println(component_io, "# trixi_commit=", commit)
    println(component_io, "# trixi_worktree=", worktree_state)
    println(component_io, "# julia_threads=", Threads.nthreads())
    println(component_io, "# command=", join(Base.ARGS, ' '))
    println(component_io, "# polydegs=", join(polydegs, ','))
    println(component_io, "# refinement_levels=", join(refinement_levels, ','))
    println(component_io, "# sigmas=", join(penalties, ','))
    println(component_io, "# auxiliary_fluxes=", join(auxiliary_fluxes, ','))
    println(component_io, "# rock4_abstol_reltol=", time_int_tol)
    println(component_io,
            "polydeg,refinement_level,ncells,h,hyperbolic_flux," *
            "auxiliary_flux,component_index,block,component," *
            "l2_error,l2_eoc,linf_error,linf_eoc")

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

            errors = run_case(polydeg, refinement_level, sigma,
                              auxiliary_flux, time_int_tol)
            error_u1 = sum(errors.l2[1:6]) / 6
            error_u2 = sum(errors.l2[7:12]) / 6
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

            for component_index in eachindex(component_names)
                block = component_index <= 6 ? "u1" : "u2"
                key = (polydeg, sigma, auxiliary_flux, component_index)
                l2_eoc = if haskey(previous_component_l2, key)
                    previous_cells, previous_error = previous_component_l2[key]
                    log(previous_error / errors.l2[component_index]) /
                    log(ncells / previous_cells)
                else
                    NaN
                end
                linf_eoc = if haskey(previous_component_linf, key)
                    previous_cells, previous_error = previous_component_linf[key]
                    log(previous_error / errors.linf[component_index]) /
                    log(ncells / previous_cells)
                else
                    NaN
                end
                format(component_io, COMPONENT_ROW_FORMAT,
                       polydeg, refinement_level, ncells, h, flux_name,
                       auxiliary_flux, component_index, block,
                       component_names[component_index],
                       errors.l2[component_index], l2_eoc,
                       errors.linf[component_index], linf_eoc)
                previous_component_l2[key] = (ncells, errors.l2[component_index])
                previous_component_linf[key] = (ncells, errors.linf[component_index])
            end
            flush(io)
            flush(component_io)
            GC.gc()
        end
    end
finally
    close(component_io)
end

println("wrote ", output_file)
println("wrote ", component_output_file)
