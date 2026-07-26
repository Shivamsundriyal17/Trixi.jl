using Plots

input_file = isempty(ARGS) ?
             joinpath(@__DIR__, "..", "results", "mms_convergence.csv") :
             abspath(ARGS[1])
output_file = length(ARGS) < 2 ?
              joinpath(@__DIR__, "..", "results", "mms_convergence.pdf") :
              abspath(ARGS[2])
mkpath(dirname(output_file))

lines = filter(line -> !isempty(line) && !startswith(line, '#'),
               readlines(input_file))
header = split(first(lines), ',')
rows = [NamedTuple{Tuple(Symbol.(header))}(Tuple(split(line, ',')))
        for line in Iterators.drop(lines, 1)]

figure = plot(layout = (1, 2), size = (1000, 430),
              xscale = :log10, yscale = :log10)
for (panel, block) in enumerate(("u1", "u2"))
    for auxiliary_flux in ("alternating", "br1"), polydeg in 1:3
        selected = filter(row -> row.block == block &&
                                     row.auxiliary_flux == auxiliary_flux &&
                                     parse(Int, row.polydeg) == polydeg, rows)
        isempty(selected) && continue
        cells = parse.(Int, getproperty.(selected, :ncells))
        errors = parse.(Float64, getproperty.(selected, :l2_error))
        label = "k=$polydeg, $auxiliary_flux"
        plot!(figure[panel], cells, errors, marker = :circle,
              linewidth = 2, label = label)
    end
    xlabel!(figure[panel], "number of cells")
    ylabel!(figure[panel], "mean L2 error")
    title!(figure[panel], block)
end

savefig(figure, output_file)
println("wrote ", output_file)
