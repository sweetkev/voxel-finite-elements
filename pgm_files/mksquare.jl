f = open("square.pgm"; write=true)

n = 4

println(f, "P2")
println(f, n, " ", n)
println(f, "255")
for j in 1:n
    for i in 1:n
        print(f, "255 ")
    end
    println(f)
end
