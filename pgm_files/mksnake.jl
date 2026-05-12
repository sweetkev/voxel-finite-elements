# This script generates a PGM file representing a snake-like pattern. The image
# is bounded by an n x n box. The snake has a specified width and gap between
# its turns. 

n = 8 # side length of bounding box
width = 2 # width of snake
vgap = 1 # vertical gap between lines of snake
output = "snake.pgm"

if length(ARGS) >= 1
    n = parse(Int, ARGS[1])
end
if length(ARGS) >= 2
    output = ARGS[2]
end

f = open(output; write=true)

nlines = div(n-width, width + vgap) + 1 # number of lines in snake

# Image header
println(f, "P2")
println(f, n, " ", n)
println(f, "255")

# Data
# First n-1 lines of snake
for line in 1:(nlines-1)
    for y in 1:(width + vgap)
        for x in 1:n
            if (y <= width)
                print(f, "255\t")
            else
                if (((line % 2 == 1) && (x > n-width)) || ((line % 2 == 0) && (x <= width)))
                    print(f,"255\t")
                else
                    print(f,"0\t")
                end
            end
        end
        println(f)
    end
end

# Final line of snake
for y in 1:width
    for x in 1:n
        print(f, "255\t")
    end
    println(f)
end

# Fill in empty space below snake
for y in 1:(n - (nlines-1)*(width + vgap) - width )
    for x in 1:n
        print(f, "0\t")
    end
    println(f)
end