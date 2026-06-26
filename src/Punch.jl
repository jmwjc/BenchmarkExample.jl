
module Punch

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 48.0
𝐷 = 12.0

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Cantilever Beam")

    Ω1,Ω2,Ω3, Γ₁, Γ₂, Γ₃, Γ₄, Γ₅ ,Γ6,Γ7,Γ8,Γ9,Γ10 = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    if transfinite > 0
        gmsh.model.mesh.setTransfiniteCurve(Γ₁, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₂, 2*transfinite-1)
        gmsh.model.mesh.setTransfiniteCurve(Γ₃, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₄,  transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₅, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ6,  2*transfinite-1)
        gmsh.model.mesh.setTransfiniteCurve(Γ7, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ8, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ9, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ10, transfinite)
        gmsh.model.mesh.setTransfiniteSurface(Ω1)
        gmsh.model.mesh.setTransfiniteSurface(Ω2)
        gmsh.model.mesh.setTransfiniteSurface(Ω3)
       
    end

    gmsh.model.mesh.setAlgorithm(2, Ω1, 1)
    gmsh.model.mesh.setAlgorithm(2, Ω2, 1)
    gmsh.model.mesh.setAlgorithm(2, Ω3, 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    gmsh.model.mesh.SecondOrderLinear = 1
    gmsh.model.mesh.secondOrderIncomplete = true
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")

    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    gmsh.model.geo.addPoint(0.0,-𝐷/2, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝐿,-𝐷/2, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  𝐿, 𝐷/2, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0, 𝐷/2, 0.0, lc, 4)

    gmsh.model.geo.addPoint(𝐷, 𝐷/2, 0.0, lc, 5)
    gmsh.model.geo.addPoint(3/4*𝐿, 𝐷/2, 0.0, lc, 6)
    gmsh.model.geo.addPoint(𝐷, -𝐷/2, 0.0, lc, 7)
    gmsh.model.geo.addPoint(3/4*𝐿, -𝐷/2, 0.0, lc, 8)

    Γ₁ = gmsh.model.geo.addLine(1, 7, 1)
    Γ₂ = gmsh.model.geo.addLine(7, 8, 2)
    Γ₃ = gmsh.model.geo.addLine(8, 2, 3)
    Γ₄ = gmsh.model.geo.addLine(2, 3, 4)
    Γ₅ = gmsh.model.geo.addLine(3, 6, 5)
    Γ6 = gmsh.model.geo.addLine(6, 5, 6)
    Γ7 = gmsh.model.geo.addLine(5, 4, 7)
    Γ8 = gmsh.model.geo.addLine(4, 1, 8)
    Γ9 = gmsh.model.geo.addLine(7, 5, 9)
    Γ10 = gmsh.model.geo.addLine(8, 6, 10)

    gmsh.model.geo.addCurveLoop([8,1,9,7],1)
    gmsh.model.geo.addCurveLoop([-9,2,10,6],2)
    gmsh.model.geo.addCurveLoop([-10,3,4,5],3)
    Ω1 = gmsh.model.geo.addPlaneSurface([1],1)
    Ω2 = gmsh.model.geo.addPlaneSurface([2],2)
    Ω3 = gmsh.model.geo.addPlaneSurface([3],3)
    gmsh.model.geo.synchronize()

    # gmsh.model.addPhysicalGroup(1, [Γ₁,Γ₃], -1, "Γʳ")
    gmsh.model.addPhysicalGroup(1, [Γ6], -1, "Γᵗ")
    gmsh.model.addPhysicalGroup(1, [Γ8,Γ₁,Γ₂,Γ₃,Γ₄ ], -1, "Γᵍ")
    gmsh.model.addPhysicalGroup(2, [Ω1,Ω2,Ω3], -1, "Ω")

    return Ω1,Ω2,Ω3, Γ₁, Γ₂, Γ₃, Γ₄, Γ₅ ,Γ6,Γ7,Γ8,Γ9,Γ10
end

end