# Run from the pinned beam environment. Results go to results/, never reference/.
using Trixi
using Printf

include(joinpath(@__DIR__, "..", "..", "test", "test_damped_intrinsic_beam_regressions.jl"))

output = isempty(ARGS) ? joinpath(@__DIR__, "results", "pre_push_validation.csv") :
         abspath(ARGS[1])
mkpath(dirname(output))

function measure_call(function_, count = 100)
    function_()
    function_()
    bytes = @allocated function_()
    seconds = minimum([@elapsed(for _ in 1:count
                                    function_()
                                end) / count for _ in 1:3])
    return bytes, seconds
end

function measure_force(force, equations)
    return measure_call(() -> force(SVector(0.2), 0.1, equations))
end

module ValidationMMS
using Trixi
end
module ValidationRotating
using Trixi
end

open(output, "w") do io
    println(io, "# julia_version=", VERSION)
    println(io, "# julia_threads=", Threads.nthreads())
    println(io, "case,ncells,l2_max,linf_max,l2_eoc,linf_eoc,seconds,allocated_bytes")
    previous = nothing
    for level in (2, 3, 4)
        seconds = @elapsed Trixi.trixi_include(ValidationMMS,
                                               joinpath(@__DIR__,
                                                        "elixir_mms_rich_physical.jl");
                                               polydeg = 2,
                                               initial_refinement_level = level,
                                               tspan = (0.0, 0.1), time_int_tol = 1.0e-12)
        errors = Base.invokelatest(ValidationMMS.analysis_callback, ValidationMMS.sol)
        l2, linf = maximum(errors.l2), maximum(errors.linf)
        l2_eoc = isnothing(previous) ? NaN : log2(previous[1] / l2)
        linf_eoc = isnothing(previous) ? NaN : log2(previous[2] / linf)
        println(io,
                join(("rich_mms", 2^level, l2, linf, l2_eoc, linf_eoc, seconds, 0), ','))
        previous = (l2, linf)
        flush(io)
    end
    bytes, seconds = Base.invokelatest(measure_force, ValidationMMS.physical_external_force,
                                       ValidationMMS.equations_hyperbolic)
    println(io, join(("forcing_call", 0, NaN, NaN, NaN, NaN, seconds, bytes), ','))
    Trixi.trixi_include(ValidationRotating, joinpath(@__DIR__, "elixir_rotating_beam.jl");
                        tspan = (0.0, 1e-4), save_times = [0.0, 1e-4])
    bytes, seconds = Base.invokelatest(measure_call,
                                       () -> Base.invokelatest(ValidationRotating.rotating_ledger_terms,
                                                               ValidationRotating.sol.u[end],
                                                               1e-4))
    println(io, join(("ledger_call", 8, NaN, NaN, NaN, NaN, seconds, bytes), ','))
end
println("Wrote ", output)
