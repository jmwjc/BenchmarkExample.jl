module SphericalShell
    
import ..BenchmarkExample
import Gmsh: gmsh

const 𝐸 = 6.825e7
const 𝜈 = 0.3
const 𝐹 = 2.0
const 𝑅 = 10.0
const 𝜃₁ = 18/180*π
const 𝜃₂ = 90/180*π
const ℎ = 0.04

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Spherical shell problem")

    Γᵇ, Γʳ, Γᵗ, Γˡ, Ω = generateGeo(lc)

    if transfinite > 0
        # transfinitex = round(transfinite*2*𝑅*𝜃/𝐿)
        # ft = floor(transfinite*2*𝑅*𝜃/𝐿)
        # fc = ceil(transfinite*2*𝑅*𝜃/𝐿)
        # transfinitex = isodd(ft) ? ft : fc
        gmsh.model.mesh.setTransfiniteCurve(Γᵇ, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γʳ, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γᵗ, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γˡ, transfinite)
        gmsh.model.mesh.setTransfiniteSurface(Ω)
    end

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    # gmsh.finalize()
end
@inline function generateGeo(lc = 1.0)
    𝑅 = 10.0
    𝜃₁ = 18/180*π # 18°
    𝜃₂ = 90/180*π # 90°
 
    𝐴 = gmsh.model.geo.addPoint(𝑅*𝜃₂, 0.0, 0.0, lc, 1)
    𝐵 = gmsh.model.geo.addPoint(0.0,  𝑅*𝜃₂, 0.0, lc, 2)
    gmsh.model.geo.addPoint(0.0, 𝑅*𝜃₁, 0.0, lc, 3)
    gmsh.model.geo.addPoint(𝑅*𝜃₁, 0.0, 0.0, lc, 4)
    Γᵇ = gmsh.model.geo.addLine(1, 2, 1)
    Γʳ = gmsh.model.geo.addLine(2, 3, 2)
    Γᵗ = gmsh.model.geo.addLine(3, 4, 3)
    Γˡ = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(0, [𝐴], -1, "𝐴")
    gmsh.model.addPhysicalGroup(0, [𝐵], -1, "𝐵")
    gmsh.model.addPhysicalGroup(1, [Γᵇ], -1, "Γᵇ")
    gmsh.model.addPhysicalGroup(1, [Γᵗ], -1, "Γᵗ")
    gmsh.model.addPhysicalGroup(1, [Γˡ], -1, "Γˡ")
    gmsh.model.addPhysicalGroup(1, [Γʳ], -1, "Γʳ")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Γᵇ, Γʳ, Γᵗ, Γˡ, Ω
end

end