module MorleysAcuteSkewPlate
    
import ..BenchmarkExample
import Gmsh: gmsh

const 𝐸 = 10.92
const 𝜈 = 0.3
const 𝐿 = 100.0
const ℎ = 0.1  # 0.1/1/10
const 𝑤 = 0.0 
const 𝜃₁ = 0.0
const 𝜃₂ = 0.0
const 𝐹 = 1.0
                   #0.1          1            10
const 𝑣ᴱ = 0.4134  #𝑣ᴱ = 0.4134  𝑣ᴱ = 0.4248  𝑣ᴱ =0.5177

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("MorleysAcuteSkewPlate")

    Γᵇ, Γʳ, Γᵗ, Γˡ, Ω = generateGeo(lc)

    if transfinite > 0
        gmsh.model.mesh.setTransfiniteCurve(Γᵇ, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γʳ, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γᵗ, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γˡ, transfinite)
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
    𝐿 = 100.0
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(𝐿, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(3^0.5*𝐿/2+𝐿, 𝐿/2 , 0.0, lc, 3)
    gmsh.model.geo.addPoint(3^0.5*𝐿/2, 𝐿/2 , 0.0, lc, 4)
    𝐴 = gmsh.model.geo.addPoint(3^0.5*𝐿/4+𝐿/2, 𝐿/4, 0.0, lc, 5)

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