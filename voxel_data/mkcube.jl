# Creates cube.nhdr (header) and cube.nrrd (raw voxel data) for a solid cube
# with n voxels along each edge. Usage: julia mkcube.jl n

n = 2
if length(ARGS) >= 1
    n = parse(Int, ARGS[1])
end

name = "cube"

# Every voxel is filled (nonzero)
data = fill(UInt8(255), n, n, n)
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
