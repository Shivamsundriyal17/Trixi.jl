module TestDampedIntrinsicBeamCampaigns

using Test, Trixi, Serialization
const beam_examples = normpath(joinpath(@__DIR__, "..", "examples",
                                        "damped_intrinsic_beam"))
const validation_output = isempty(ARGS) ? mktempdir() : abspath(first(ARGS))
mkpath(validation_output)
function run_script(name, arguments; env = ())
    previous_args = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, arguments)
    try
        return withenv(env...) do
            sandbox = Module(gensym(:BeamCampaign))
            Base.include(sandbox, joinpath(beam_examples, name))
            sandbox
        end
    finally
        empty!(ARGS)
        append!(ARGS, previous_args)
    end
end
@testset "Rotating campaign generation and reuse" begin
    folder = joinpath(validation_output, "rotating")
    run_script("run_rotating_beam_campaign.jl", [folder];
               env = ("ROTATING_CAMPAIGN_QUICK" => "true",
                      "ROTATING_REUSE_RESULTS" => "false"))
    run_script("run_rotating_beam_campaign.jl", [folder];
               env = ("ROTATING_CAMPAIGN_QUICK" => "true",
                      "ROTATING_REUSE_RESULTS" => "true"))
    caught = try
        run_script("run_rotating_beam_campaign.jl", [folder];
                   env = ("ROTATING_CAMPAIGN_QUICK" => "false",
                          "ROTATING_REUSE_RESULTS" => "true"))
        nothing
    catch error
        error
    end
    @test !isnothing(caught)
    @test occursin("configuration differs", sprint(showerror, caught))
    @test isfile(joinpath(folder, "summary.csv"))
end
@testset "Cantilever checkpoint resume" begin
    output = joinpath(validation_output, "zero_load_sweep.csv")
    environment = ("CANTILEVER_ACCELERATION_RMS_G" => "0.0",
                   "CANTILEVER_GRAVITY_MULTIPLIER" => "0.0",
                   "CANTILEVER_POLYDEG" => "2", "CANTILEVER_REFINEMENT_LEVEL" => "1",
                   "CANTILEVER_SAMPLES_PER_CYCLE" => "8",
                   "CANTILEVER_SWEEP_NORMALIZED_FREQUENCIES" => "1.0",
                   "CANTILEVER_SWEEP_FIRST_CYCLES" => "2",
                   "CANTILEVER_SWEEP_CONTINUATION_CYCLES" => "2",
                   "CANTILEVER_SWEEP_ADDITIONAL_CYCLES" => "2",
                   "CANTILEVER_SWEEP_MAXIMUM_CYCLES" => "2",
                   "CANTILEVER_SWEEP_RAMP_CYCLES" => "0",
                   "CANTILEVER_SWEEP_RELTOL" => "1e-4",
                   "CANTILEVER_SWEEP_OUTPUT" => output)
    run_script("run_base_excited_cantilever_sweep.jl", String[];
               env = (environment..., "CANTILEVER_SWEEP_RESUME" => "false"))
    original = read(output)
    run_script("run_base_excited_cantilever_sweep.jl", String[];
               env = (environment..., "CANTILEVER_SWEEP_RESUME" => "true"))
    @test read(output) == original
    checkpoint = deserialize(replace(output, ".csv" => "_checkpoint.jls"))
    @test checkpoint.format_version == 5
    caught = try
        changed = filter(pair -> first(pair) != "CANTILEVER_SWEEP_RELTOL", environment)
        run_script("run_base_excited_cantilever_sweep.jl", String[];
                   env = (changed..., "CANTILEVER_SWEEP_RESUME" => "true",
                          "CANTILEVER_SWEEP_RELTOL" => "1e-3"))
        nothing
    catch error
        error
    end
    @test !isnothing(caught)
    @test occursin("configuration differs for reltol", sprint(showerror, caught))
    @test read(output) == original
end

@testset "Zero-gravity checkpoint cycle replay" begin
    # Construct a portable, zero-gravity equilibrium checkpoint in legacy format.
    setup = run_script("elixir_base_excited_cantilever.jl", String[];
                       env = ("CANTILEVER_SETUP_ONLY" => "true",
                              "CANTILEVER_GRAVITY_MULTIPLIER" => "0.0",
                              "CANTILEVER_ACCELERATION_RMS_G" => "0.0",
                              "CANTILEVER_POLYDEG" => "2",
                              "CANTILEVER_REFINEMENT_LEVEL" => "1"))
    checkpoint = joinpath(validation_output, "zero_gravity_checkpoint.jls")
    row = (; frequency_hz = setup.frequency_hz, omega = setup.omega,
           normalized_frequency = 1.0)
    serialize(checkpoint,
              (; format_version = 4, acceleration_rms_g = 0.0,
               gravity_multiplier = 0.0, polydeg = setup.polydeg,
               refinement_level = setup.refinement_level,
               constraint_multiplier = setup.constraint_multiplier,
               rows = [row], continuation_state = copy(setup.split_problem.u0)))
    sandbox = run_script("save_base_excited_cantilever_cycle.jl",
                         [checkpoint, joinpath(validation_output, "zero_gravity_cycle")];
                         env = ("CANTILEVER_CYCLE_FRAMES" => "8",
                                "CANTILEVER_GEOMETRY_POINTS_PER_ELEMENT" => "8",
                                "CANTILEVER_GRAVITY_MULTIPLIER" => "1.0"))
    @test sandbox.gravity_multiplier == 0.0
    @test sandbox.rhs_parameters.gamma == 0.0
    @test sandbox.cycle_solution.t[end] == sandbox.cycle_period
end
println("CAMPAIGN_INTEGRATION_COMPLETE")

end # module
