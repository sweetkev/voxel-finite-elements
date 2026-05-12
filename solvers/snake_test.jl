#!/usr/bin/env julia

using Printf

m_values = isempty(ARGS) ? collect(3:12) : [parse(Int, x) for x in ARGS]
nlevels = 4
src_pgm = joinpath("..", "pgm_files", "snake.pgm")
mksnake = joinpath("..", "pgm_files", "mksnake.jl")
dirichlet_exe = joinpath(".", "dirichlet_2d")
graph_exe = joinpath(".", "graph_solver")

function run_and_capture(cmd::Cmd)
    println("Running: ", cmd)
    return read(cmd, String)
end

function parse_cg_iterations(output::AbstractString)
    patterns = [r"(?i)(?:cg|conjugate gradient).*?(\d+)", r"(\d+)\s*iterations", r"iterations.*?(\d+)", r"(?i)iteration\s*:\s*(\d+)"]
    all_iters = Int[]
    for pat in patterns
        for m in eachmatch(pat, output)
            push!(all_iters, parse(Int, m.captures[1]))
        end
    end
    if isempty(all_iters)
        return missing
    else
        return maximum(all_iters) + 1  # since iterations start at 0
    end
end

function ensure_executable(path::String, source::String)
    if !isfile(path)
        println("Executable $(path) not found. Attempting to compile $(source).")
        run(`g++ -O3 $(source) -o $(path)`)    
    end
end

ensure_executable(dirichlet_exe, "dirichlet_2d.cpp")
ensure_executable(graph_exe, "graph_solver.cpp")

f = open("snake_test.tex"; write=true)
println(f,"\\documentclass{article}")
println(f,"\\begin{document}")
println(f,"Number of levels: $nlevels\\\\")
println(f,"\\begin{tabular}{ccc}")
println(f,"\\hline")
println(f,"n & na\\\"ive CG iter & graph CG iter \\\\")
println(f,"\\hline")
for m in m_values
    n = Int(1) << m
    meshfile = src_pgm
    println("\nGenerating mesh for m=$m, n=2^m=$n -> $meshfile")
    run(`julia $(mksnake) $(n) $(meshfile)`)

    dir_output = run_and_capture(`$(dirichlet_exe) -f $(meshfile) --nlevels $(nlevels)`)
    graph_output = run_and_capture(`$(graph_exe) -f $(meshfile) --nlevels $(nlevels)`)

    dir_iter = parse_cg_iterations(dir_output)
    graph_iter = parse_cg_iterations(graph_output)

    dir_str = dir_iter === missing ? "?" : string(dir_iter)
    graph_str = graph_iter === missing ? "?" : string(graph_iter)

    println(f,"\$2^{$(m)}\$ & $(dir_str) & $(graph_str) \\\\")
end
println(f,"\\hline")
println(f,"\\end{tabular}")
println(f,"\\end{document}")
