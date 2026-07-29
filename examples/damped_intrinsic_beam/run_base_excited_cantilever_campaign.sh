#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 {02g_upper|02g_lower|05g_upper|05g_lower|02g_refined_lower|05g_refined_upper}" >&2
    exit 2
fi

example_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repository_directory="$(cd -- "${example_directory}/../.." && pwd)"
results_directory="${example_directory}/results"
mkdir -p "${results_directory}"

acceleration=
polydeg=4
refinement_level=1
frequencies=
first_cycles=
continuation_cycles=
additional_cycles=
maximum_cycles=
periodicity_tolerance=
output_stem=

case "$1" in
    02g_upper)
        acceleration=0.2
        frequencies="0.8231814788683266,0.848887538640635,0.8745935845027772,0.9003017515931105,0.9260088580693485,0.9517159539394453,0.9774220104077285,0.9928458500478329,1.0031291062778251,1.0082707429999171,1.0103273774149031,1.012382975732101,1.0144396101470867,1.0164962690783801,1.0185529140995069,1.0206085018105635,1.0309244496373067,1.0330190508991888,1.0351125564245272,1.0372060864661727,1.0393028161192599,1.041395274940669,1.0424425596623985,1.047678987269137,1.0529133320741579,1.0633862646955452,1.073856032827939,1.1000341144845291,1.1262132323780092,1.1523902570027325"
        first_cycles=40
        continuation_cycles=18
        additional_cycles=10
        maximum_cycles=70
        periodicity_tolerance=0.002
        output_stem=farokhi_02g_k4n2_sweep
        ;;
    02g_lower)
        acceleration=0.2
        frequencies="1.1523902570027325,1.1262132323780092,1.1000341144845291,1.073856032827939,1.0633862646955452,1.0529133320741579,1.047678987269137,1.0424425596623985,1.041395274940669,1.0393028161192599,1.0372060864661727,1.0351125564245272,1.0330190508991888,1.0309244496373067"
        first_cycles=80
        continuation_cycles=30
        additional_cycles=20
        maximum_cycles=190
        periodicity_tolerance=0.0002
        output_stem=farokhi_02g_k4n2_down_sweep
        ;;
    05g_upper)
        acceleration=0.5
        frequencies="0.8231814788683266,0.848887538640635,0.8745935845027772,0.9003017515931105,0.9260088580693485,0.9517159539394453,0.9774220104077285,1.0031291062778251,1.0185529459179299,1.0236935432382079,1.0288351627461085,1.0308928822914973,1.0319201511430482,1.0329495307557512,1.0350051256298216,1.0370617740940755,1.0391205403433936,1.0411750674403546,1.0432306552905128"
        first_cycles=80
        continuation_cycles=30
        additional_cycles=20
        maximum_cycles=190
        periodicity_tolerance=0.0005
        output_stem=farokhi_05g_k4n2_up_sweep
        ;;
    05g_lower)
        acceleration=0.5
        frequencies="1.131663486450082,1.1059574268168753,1.0802492807997222,1.059683821651717,1.0545421849296253,1.0494016121256542,1.0473449671045272,1.0452872614693045,1.0432306130050508,1.0411750146878531,1.0391205086640722,1.037061721341574,1.0350050728773201,1.0329494990764299"
        first_cycles=80
        continuation_cycles=30
        additional_cycles=20
        maximum_cycles=190
        periodicity_tolerance=0.0002
        output_stem=farokhi_05g_k4n2_down_sweep
        ;;
    02g_refined_lower)
        acceleration=0.2
        polydeg=3
        refinement_level=2
        frequencies="1.041395274940669"
        first_cycles=180
        continuation_cycles=180
        additional_cycles=20
        maximum_cycles=260
        periodicity_tolerance=0.0002
        output_stem=farokhi_02g_k3n4_refined_check
        ;;
    05g_refined_upper)
        acceleration=0.5
        polydeg=3
        refinement_level=2
        frequencies="0.9774220104077285,1.0031291062778251,1.0185529459179299,1.0288351627461085,1.0350051256298216,1.0391205403433936,1.0411750674403546"
        first_cycles=130
        continuation_cycles=20
        additional_cycles=20
        maximum_cycles=190
        periodicity_tolerance=0.0005
        output_stem=farokhi_05g_k3n4_refined_up_sweep
        ;;
    *)
        echo "unknown campaign: $1" >&2
        exit 2
        ;;
esac

env \
    CANTILEVER_ACCELERATION_RMS_G="${acceleration}" \
    CANTILEVER_POLYDEG="${polydeg}" \
    CANTILEVER_REFINEMENT_LEVEL="${refinement_level}" \
    CANTILEVER_SAMPLES_PER_CYCLE="${CANTILEVER_SAMPLES_PER_CYCLE:-96}" \
    CANTILEVER_SWEEP_NORMALIZED_FREQUENCIES="${frequencies}" \
    CANTILEVER_SWEEP_FIRST_CYCLES="${first_cycles}" \
    CANTILEVER_SWEEP_CONTINUATION_CYCLES="${continuation_cycles}" \
    CANTILEVER_SWEEP_ADDITIONAL_CYCLES="${additional_cycles}" \
    CANTILEVER_SWEEP_MAXIMUM_CYCLES="${maximum_cycles}" \
    CANTILEVER_SWEEP_PERIODICITY_TOLERANCE="${periodicity_tolerance}" \
    CANTILEVER_SWEEP_ARCHIVE_STATES="${CANTILEVER_SWEEP_ARCHIVE_STATES:-true}" \
    CANTILEVER_SWEEP_RESUME="${CANTILEVER_SWEEP_RESUME:-true}" \
    CANTILEVER_SWEEP_OUTPUT="${results_directory}/${output_stem}.csv" \
    CANTILEVER_SWEEP_CHECKPOINT="${results_directory}/${output_stem}_checkpoint.jls" \
    julia --compiled-modules=no \
        --project="${example_directory}" \
        "${example_directory}/run_base_excited_cantilever_sweep.jl"
