# Longer local PR checks; run with the pinned beam environment and one Julia thread.
using Trixi
using Test

output = joinpath(@__DIR__, "results", "extended_validation.txt")
mkpath(dirname(output))
open(output, "w") do io
    println(io, "Julia ", VERSION, "; threads=", Threads.nthreads())
    for steady in (false, true)
        sandbox = Module(gensym(:ExtendedRotating))
        elapsed = @elapsed withenv("ROTATING_T_END" => "2.0",
                                   "ROTATING_STEADY_INITIAL" => string(steady),
                                   "ROTATING_POLYDEG" => "3",
                                   "ROTATING_REFINEMENT_LEVEL" => "3",
                                   "ROTATING_SAVE_COUNT" => "101",
                                   "ROTATING_RECORD_LEDGER" => "true",
                                   "ROTATING_OUTPUT_FILE" => joinpath(@__DIR__, "results",
                                                                     "extended_rotating_$(steady).jls")) do
            Base.include(sandbox, joinpath(@__DIR__, "save_rotating_beam_data.jl"))
        end
        @test sandbox.sol.t[end] == 2.0
        @test sandbox.maximum_ledger_relative_residual < 1.0e-10
        println(io, "rotating steady=", steady, "; seconds=", elapsed,
                "; steps=", sandbox.sol.destats.naccept,
                "; diagnostics=", sandbox.diagnostics)
        flush(io)
        println("Completed rotating steady=", steady, " in ", elapsed, " seconds")
        flush(stdout)
    end
    sandbox = Module(gensym(:ExtendedCantilever))
    elapsed = @elapsed withenv("CANTILEVER_SETUP_ONLY" => "false",
                               "CANTILEVER_ACCELERATION_RMS_G" => "0.2",
                               "CANTILEVER_GRAVITY_MULTIPLIER" => "1.0",
                               "CANTILEVER_FREQUENCY_HZ" => "9.2",
                               "CANTILEVER_POLYDEG" => "4",
                               "CANTILEVER_REFINEMENT_LEVEL" => "1",
                               "CANTILEVER_RAMP_CYCLES" => "8",
                               "CANTILEVER_SETTLING_CYCLES" => "10",
                               "CANTILEVER_MEASUREMENT_CYCLES" => "2",
                               "CANTILEVER_SAMPLES_PER_CYCLE" => "96") do
        Base.include(sandbox, joinpath(@__DIR__, "elixir_base_excited_cantilever.jl"))
    end
    @test sandbox.solution.t[end] == sandbox.t_end
    @test sandbox.result.transverse_peak > 0
    ledger = Base.invokelatest(sandbox.cantilever_cycle_ledger,
                              sandbox.solution.u, sandbox.solution.t)
    @test isfinite(ledger.relative_ledger_residual)
    println(io, "cantilever cycles=20; seconds=", elapsed,
            "; result=", sandbox.result, "; ledger=", ledger)
    println("Completed driven cantilever in ", elapsed, " seconds")
end
println("EXTENDED_VALIDATION_COMPLETE: ", output)
