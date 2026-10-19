# Creates spacefilling.nhdr (header) and spacefilling.raw (raw voxel data) for
# a 3D Hilbert curve of a given order. The curve visits a 2^order x 2^order x
# 2^order grid of points, each drawn as a cube of voxels with a specified width
# and gap between neighboring points. The resulting box has side length
# n = 2^order*(width + gap) - gap. Usage: julia mkspacefilling.jl [order] [name]

order = 2 # order of Hilbert curve
width = 2 # width of curve
gap = 1 # gap between adjacent segments of curve
name = "spacefilling"

if length(ARGS) >= 1
    order = parse(Int, ARGS[1])
end
if length(ARGS) >= 2
    name = ARGS[2]
end

m = 2^order # number of curve points along each edge
n = m*(width + gap) - gap # side length of bounding box

# Returns the 0-based grid coordinates of the h-th point (0-based) on the
# Hilbert curve, using Skilling's transpose algorithm.
function hilbert_point(h, order)
    # Distribute bits of h into transposed form, most significant bit first
    X = zeros(Int, 3)
    for b in 0:(3*order - 1)
        bit = (h >> (3*order - 1 - b)) & 1
        X[b % 3 + 1] |= bit << (order - 1 - div(b, 3))
    end

    # Gray decode
    t = X[3] >> 1
    for i in 3:-1:2
        X[i] ⊻= X[i-1]
    end
    X[1] ⊻= t

    # Undo excess work
    Q = 2
    while Q != (1 << order)
        P = Q - 1
        for i in 3:-1:1
            if X[i] & Q != 0
                X[1] ⊻= P
            else
                t = (X[1] ⊻ X[i]) & P
                X[1] ⊻= t
                X[i] ⊻= t
            end
        end
        Q <<= 1
    end
    return X
end

# Voxel index range covering grid points from a to b along one axis
voxels(a, b) = (min(a, b)*(width + gap) + 1):(max(a, b)*(width + gap) + width)

# Draw each point and the segment connecting it to the previous point
data = zeros(UInt8, n, n, n)
prev = hilbert_point(0, order)
for h in 0:(m^3 - 1)
    global prev
    p = hilbert_point(h, order)
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
