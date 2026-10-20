# Creates spacefilling.nhdr (header) and spacefilling.raw (raw voxel data) for
# a space-filling curve built from stacked snakes. The curve visits an m x m x m
# grid of points, each drawn as a cube of voxels with a specified width and gap
# between neighboring points. On each z-level the curve traces the snake pattern
# from pgm_files/mksnake.jl, then steps up one level and traces the snake in
# reverse, so it ends directly above where the previous level ended. The
# resulting box has side length n = m*(width + gap) - gap.
# Usage: julia mkspacefilling.jl [m] [name]

m = 4 # number of "curve" points along each edge
width = 3 # width of curve
gap = 1 # gap between adjacent segments of curve
name = "spacefilling"

if length(ARGS) >= 1
    m = parse(Int, ARGS[1])
end
if length(ARGS) >= 2
    name = ARGS[2]
end

n = m*(width + gap) - gap # side length of bounding box

# 0-based (x, y) grid points of the 2D snake: each row runs along x, alternating
# direction, and rows advance along y
snake = [(isodd(y) ? m - 1 - x : x, y) for y in 0:(m - 1) for x in 0:(m - 1)]

# Stack the snake on each level, reversing it on odd levels so each level
# starts directly above where the previous one ended
points = [[p[1], p[2], z] for z in 0:(m - 1) for p in (iseven(z) ? snake : reverse(snake))]

# Voxel index range covering grid points from a to b along one axis
voxels(a, b) = (min(a, b)*(width + gap) + 1):(max(a, b)*(width + gap) + width)

# Draw each point and the segment connecting it to the previous point
data = zeros(UInt8, n, n, n)
prev = points[1]
for p in points
    global prev
    data[voxels(prev[1], p[1]), voxels(prev[2], p[2]), voxels(prev[3], p[3])] .= 255
    prev = p
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
