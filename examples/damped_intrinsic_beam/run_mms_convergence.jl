using Dates: now, UTC
using Printf: Format, @printf, @sprintf, format
using Trixi

const ELIXIR_OVERRIDE = get(ENV, "MMS_ELIXIR", "")
const ELIXIR = isempty(ELIXIR_OVERRIDE) ?
               joinpath(@__DIR__, "elixir_mms_rich_physical.jl") :
               abspath(ELIXIR_OVERRIDE)
const IS_RICH_MMS = basename(ELIXIR) == "elixir_mms_rich_physical.jl"
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const COMPONENT_ROW_FORMAT = Format("%d,%d,%d,%.16e,%s,%s,%d,%s,%s," *
                                    "%.16e,%.8f,%.16e,%.8f\n")
const VISCOUS_ROW_FORMAT = Format("%d,%d,%d,%.16e,%s,%s,%d,%s," *
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
viscous_output_file = length(ARGS) < 3 ?
                      first(splitext(output_file)) *
                      "_viscous_componentwise.csv" :
                      abspath(ARGS[3])
component_names = ("v1", "v2", "v3", "omega1", "omega2", "omega3",
                   "f1", "f2", "f3", "m1", "m2", "m3")
viscous_component_names = ("r_tau_f1", "r_tau_f2", "r_tau_f3",
                           "r_tau_m1", "r_tau_m2", "r_tau_m3")

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

function viscous_resultant_error_norms(semi, sol, analyzer,
                                       manufactured_solution_t)
    final_state = sol.u[end]
    final_time = sol.t[end]

    # Re-evaluate the actual parabolic operator at the final numerical state.
    # This populates the LDG/BR1 auxiliary gradient stored by Trixi without
    # altering the state used by the evolution.
    du = similar(final_state)
    Base.invokelatest(Trixi.rhs_parabolic!, du, final_state, semi, final_time)

    mesh = semi.mesh
    equations_parabolic = semi.equations_parabolic
    equations_hyperbolic = equations_parabolic.equations_hyperbolic
    solver = semi.solver
    cache = semi.cache
    cache_parabolic = semi.cache_parabolic
    u = Trixi.wrap_array(final_state, mesh, equations_parabolic, solver,
                         cache_parabolic)
    gradients = cache_parabolic.viscous_container.gradients

    vandermonde = analyzer.vandermonde
    weights = analyzer.weights
    node_coordinates = cache.elements.node_coordinates
    n_analysis_nodes = Trixi.nnodes(analyzer)
    u_local = Array{eltype(final_state)}(undef, 12, n_analysis_nodes)
    gradient_local = similar(u_local)
    x_local = Array{eltype(node_coordinates)}(undef, 1, n_analysis_nodes)
    l2_error = zeros(eltype(final_state), 6)
    linf_error = zeros(eltype(final_state), 6)

    for element in Trixi.eachelement(solver, cache)
        Trixi.multiply_dimensionwise!(u_local, vandermonde,
                                      view(u, :, :, element))
        Trixi.multiply_dimensionwise!(gradient_local, vandermonde,
                                      view(gradients, :, :, element))
        Trixi.multiply_dimensionwise!(x_local, vandermonde,
                                      view(node_coordinates, :, :, element))
        volume_jacobian = Trixi.volume_jacobian(element, mesh, cache)

        for i in Trixi.eachnode(analyzer)
            u_node = SVector{12}(ntuple(variable -> u_local[variable, i], 12))
            gradient_node = SVector{12}(ntuple(variable ->
                                               gradient_local[variable, i], 12))
            r_tau_numerical = intrinsic_beam_damping_resultant(u_node,
                                                               gradient_node,
                                                               equations_parabolic)
            exact_time_derivative = Base.invokelatest(manufactured_solution_t,
                                                       x_local[1, i], final_time)
            exact_u2_t = SVector{6}(ntuple(component ->
                                           exact_time_derivative[6 + component], 6))
            r_tau_exact = equations_hyperbolic.damping_matrix * exact_u2_t
            difference = r_tau_exact - r_tau_numerical
            quadrature_weight = weights[i] * volume_jacobian
            l2_error .+= difference .^ 2 .* quadrature_weight
            linf_error .= max.(linf_error, abs.(difference))
        end
    end

    l2_error .= sqrt.(l2_error ./ Trixi.total_volume(mesh))
    return (l2 = l2_error, linf = linf_error)
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
    if get(ENV, "MMS_PRINT_EXACT_AUDIT", "0") == "1" &&
       isdefined(MMSRun, :rich_mms_exact_audit)
        println("rich_mms_exact_audit=", MMSRun.rich_mms_exact_audit)
        println("rich_mms_boundary_audit=", MMSRun.rich_mms_boundary_audit)
    end
    errors = Base.invokelatest(MMSRun.analysis_callback, MMSRun.sol)
    viscous_errors = viscous_resultant_error_norms(MMSRun.semi, MMSRun.sol,
                                                   MMSRun.analysis_callback.affect!.analyzer,
                                                   MMSRun.manufactured_solution_t)
    return (l2 = collect(errors.l2), linf = collect(errors.linf),
            viscous_l2 = collect(viscous_errors.l2),
            viscous_linf = collect(viscous_errors.linf))
end

mkpath(dirname(output_file))
mkpath(dirname(component_output_file))
mkpath(dirname(viscous_output_file))
commit, worktree_state = repository_metadata()
previous = Dict{Tuple{Int, Float64, String, String}, Tuple{Int, Float64}}()
previous_component_l2 = Dict{Tuple{Int, Float64, String, Int},
                             Tuple{Int, Float64}}()
previous_component_linf = Dict{Tuple{Int, Float64, String, Int},
                               Tuple{Int, Float64}}()
previous_viscous_l2 = Dict{Tuple{Int, Float64, String, Int}, Tuple{Int, Float64}}()
previous_viscous_linf = Dict{Tuple{Int, Float64, String, Int}, Tuple{Int, Float64}}()

component_io = open(component_output_file, "w")
viscous_io = open(viscous_output_file, "w")
try
    println(component_io,
            "# Componentwise physical Kelvin-Voigt MMS errors")
    println(component_io, "# generated_utc=", now(UTC))
    println(component_io, "# julia_version=", VERSION)
    println(component_io, "# trixi_commit=", commit)
    println(component_io, "# trixi_worktree=", worktree_state)
    println(component_io, "# julia_threads=", Threads.nthreads())
    println(component_io, "# command=", join(Base.ARGS, ' '))
    println(component_io, "# elixir=", ELIXIR)
    IS_RICH_MMS && println(component_io, "# rich_mms_lambda=",
                           get(ENV, "RICH_MMS_LAMBDA", "2.0"))
    IS_RICH_MMS && println(component_io, "# rich_mms_profile=",
                           get(ENV, "RICH_MMS_PROFILE", "exp1"))
    println(component_io, "# polydegs=", join(polydegs, ','))
    println(component_io, "# refinement_levels=", join(refinement_levels, ','))
    println(component_io, "# sigmas=", join(penalties, ','))
    println(component_io, "# auxiliary_fluxes=", join(auxiliary_fluxes, ','))
    println(component_io, "# rock4_abstol_reltol=", time_int_tol)
    println(component_io,
            "polydeg,refinement_level,ncells,h,hyperbolic_flux," *
            "auxiliary_flux,component_index,block,component," *
            "l2_error,l2_eoc,linf_error,linf_eoc")

    println(viscous_io, "# Componentwise numerical Kelvin-Voigt resultant errors")
    println(viscous_io, "# r_tau is reconstructed from Trixi's discrete auxiliary gradient")
    println(viscous_io, "# generated_utc=", now(UTC))
    println(viscous_io, "# julia_version=", VERSION)
    println(viscous_io, "# trixi_commit=", commit)
    println(viscous_io, "# trixi_worktree=", worktree_state)
    println(viscous_io, "# julia_threads=", Threads.nthreads())
    println(viscous_io, "# command=", join(Base.ARGS, ' '))
    println(viscous_io, "# elixir=", ELIXIR)
    IS_RICH_MMS && println(viscous_io, "# rich_mms_lambda=",
                           get(ENV, "RICH_MMS_LAMBDA", "2.0"))
    IS_RICH_MMS && println(viscous_io, "# rich_mms_profile=",
                           get(ENV, "RICH_MMS_PROFILE", "exp1"))
    println(viscous_io, "# polydegs=", join(polydegs, ','))
    println(viscous_io, "# refinement_levels=", join(refinement_levels, ','))
    println(viscous_io, "# sigmas=", join(penalties, ','))
    println(viscous_io, "# auxiliary_fluxes=", join(auxiliary_fluxes, ','))
    println(viscous_io, "# rock4_abstol_reltol=", time_int_tol)
    println(viscous_io,
            "polydeg,refinement_level,ncells,h,hyperbolic_flux," *
            "auxiliary_flux,component_index,component," *
            "l2_error,l2_eoc,linf_error,linf_eoc")

    open(output_file, "w") do io
        println(io, "# Physical, constitutively consistent Kelvin-Voigt MMS")
        println(io, "# generated_utc=", now(UTC))
        println(io, "# julia_version=", VERSION)
        println(io, "# trixi_commit=", commit)
        println(io, "# trixi_worktree=", worktree_state)
        println(io, "# julia_threads=", Threads.nthreads())
        println(io, "# command=", join(Base.ARGS, ' '))
        println(io, "# elixir=", ELIXIR)
        IS_RICH_MMS && println(io, "# rich_mms_lambda=",
                               get(ENV, "RICH_MMS_LAMBDA", "2.0"))
        IS_RICH_MMS && println(io, "# rich_mms_profile=",
                               get(ENV, "RICH_MMS_PROFILE", "exp1"))
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
            error_r_tau = sum(errors.viscous_l2) / 6
            @printf("u1=%.8e, u2=%.8e, r_tau=%.8e\n", error_u1, error_u2,
                    error_r_tau)

            for (block, error) in (("u1", error_u1), ("u2", error_u2),
                                   ("r_tau", error_r_tau))
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
            for component_index in eachindex(viscous_component_names)
                key = (polydeg, sigma, auxiliary_flux, component_index)
                l2_eoc = if haskey(previous_viscous_l2, key)
                    previous_cells, previous_error = previous_viscous_l2[key]
                    log(previous_error / errors.viscous_l2[component_index]) /
                    log(ncells / previous_cells)
                else
                    NaN
                end
                linf_eoc = if haskey(previous_viscous_linf, key)
                    previous_cells, previous_error = previous_viscous_linf[key]
                    log(previous_error / errors.viscous_linf[component_index]) /
                    log(ncells / previous_cells)
                else
                    NaN
                end
                format(viscous_io, VISCOUS_ROW_FORMAT,
                       polydeg, refinement_level, ncells, h, flux_name,
                       auxiliary_flux, component_index,
                       viscous_component_names[component_index],
                       errors.viscous_l2[component_index], l2_eoc,
                       errors.viscous_linf[component_index], linf_eoc)
                previous_viscous_l2[key] = (ncells,
                                            errors.viscous_l2[component_index])
                previous_viscous_linf[key] = (ncells,
                                              errors.viscous_linf[component_index])
            end
            flush(io)
            flush(component_io)
            flush(viscous_io)
            GC.gc()
        end
    end
finally
    close(component_io)
    close(viscous_io)
end

println("wrote ", output_file)
println("wrote ", component_output_file)
println("wrote ", viscous_output_file)
