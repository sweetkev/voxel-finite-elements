f = open("square.pgm"; write=true)

n = 2
if length(ARGS) >= 1
    n = parse(Int, ARGS[1])
end

println(f, "P2")
println(f, n, " ", n)
println(f, "255")
for j in 1:n
    for i in 1:n
        print(f, "255 ")
    end
    println(f)
end
