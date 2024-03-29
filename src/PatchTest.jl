
module PatchTest

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 1.0

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Patch Test")

    Ω, Γ₁, Γ₂, Γ₃, Γ₄ = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    if transfinite > 0
        gmsh.model.mesh.setTransfiniteCurve(Γ₁, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₂, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₃, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₄, transfinite)
        gmsh.model.mesh.setTransfiniteSurface(Ω)
    end

    gmsh.model.mesh.setAlgorithm(2, Ω, 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝐿, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  𝐿,   𝐿, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0,   𝐿, 0.0, lc, 4)
    Γ₁ = gmsh.model.geo.addLine(1, 2, 1)
    Γ₂ = gmsh.model.geo.addLine(2, 3, 2)
    Γ₃ = gmsh.model.geo.addLine(3, 4, 3)
    Γ₄ = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [Γ₁], -1, "Γ¹")
    gmsh.model.addPhysicalGroup(1, [Γ₂], -1, "Γ²")
    gmsh.model.addPhysicalGroup(1, [Γ₃], -1, "Γ³")
    gmsh.model.addPhysicalGroup(1, [Γ₄], -1, "Γ⁴")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Ω, Γ₁, Γ₂, Γ₃, Γ₄
end
end