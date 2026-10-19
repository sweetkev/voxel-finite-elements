# Creates folds.nhdr (header) and folds.raw (raw voxel data) for a snake-like
# pattern of folds bounded by an n x n x n box. Each slice in z matches the
# image produced by pgm_files/mksnake.jl. The folds have a specified width and
# gap between turns. Usage: julia mkfolds.jl [n] [name]

n = 8 # side length of bounding box
width = 2 # width of folds
vgap = 1 # gap between folds
name = "folds"

if length(ARGS) >= 1
    n = parse(Int, ARGS[1])
end
if length(ARGS) >= 2
    name = ARGS[2]
end

nlines = div(n-width, width + vgap) + 1 # number of folds

# Returns true if pixel (x, y) of the 2D snake pattern is filled
function filled(x, y)
    # Empty space below the final fold
    if y > (nlines-1)*(width + vgap) + width
        return false
    end
    line = div(y-1, width + vgap) + 1
    offset = (y-1) % (width + vgap) + 1
    if offset <= width
        return true
    end
    # Connector between folds, alternating sides
    return ((line % 2 == 1) && (x > n-width)) || ((line % 2 == 0) && (x <= width))
end

# Extrude the 2D pattern along z (x varies fastest, then y, then z)
data = zeros(UInt8, n, n, n)
for z in 1:n, y in 1:n, x in 1:n
    if filled(x, y)
        data[x, y, z] = 255
    end
end

open(name * ".raw"; write=true) do f
    write(f, data)
end

open(name * ".nhdr"; write=true) do f
    println(f, "NRRD0004")
    println(f, "# Complete NRRD file format specification at:")
    println(f, "# http://teem.sourceforge.net/nrrd/format.html")
    println(f, "type: uchar")
    println(f, "endian: little")
    println(f, "encoding: raw")
    println(f, "dimension: 3")
    println(f, "space: 3D-right-handed")
    println(f, "kinds: domain domain domain")
    println(f, "sizes: ", n, " ", n, " ", n)
    println(f, "space directions: (1,0,0) (0,1,0) (0,0,1)")
    println(f, "space origin: (0.0,0.0,0.0)")
    println(f, "data file: ", name, ".raw")
end
