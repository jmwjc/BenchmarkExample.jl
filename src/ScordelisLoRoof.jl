module ScordelisLoRoof
    
import ..BenchmarkExample
import Gmsh: gmsh

const 𝐸 = 4.32e8
const 𝜈 = 0.0
const 𝑏₂ = 90.0
const 𝑣ₐ = 0.3024
const 𝑣ₘ = 0.30078086

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Scordelis-Lo roof problem")

    Γᵇ, Γʳ, Γᵗ, Γˡ, Ω = generateGeo(lc)

    if transfinite > 0
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
    dim, tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(dim, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    𝑅 = 25.0
    𝐿 = 50.0
    𝜃 = 40/180*π # 40°
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(𝑅*𝜃, 0.0, 0.0, lc, 2)
    𝐴 = gmsh.model.geo.addPoint(𝑅*𝜃,  𝐿/2, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0, 𝐿/2, 0.0, lc, 4)
    Γᵇ = gmsh.model.geo.addLine(1, 2, 1)
    Γʳ = gmsh.model.geo.addLine(2, 3, 2)
    Γᵗ = gmsh.model.geo.addLine(3, 4, 3)
    Γˡ = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(0, [𝐴], -1, "𝐴")
    gmsh.model.addPhysicalGroup(1, [Γᵇ], -1, "Γᵇ")
    gmsh.model.addPhysicalGroup(1, [Γᵗ], -1, "Γᵗ")
    gmsh.model.addPhysicalGroup(1, [Γˡ], -1, "Γˡ")
    gmsh.model.addPhysicalGroup(1, [Γʳ], -1, "Γʳ")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Γᵇ, Γʳ, Γᵗ, Γˡ, Ω
end

end