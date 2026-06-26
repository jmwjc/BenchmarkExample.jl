
module CookMembrance

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 48.0
𝐷 = 44
d = 16

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false,coef = 1.0)
    gmsh.initialize()
    gmsh.model.add("Cook membrance")

    Ω, Γ₁, Γ₂, Γ₃, Γ₄ = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    if transfinite > 0
        # gmsh.model.mesh.setTransfiniteCurve(Γ₁, 2*transfinite-1, "Progression",-coef)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₁, 4*transfinite-3)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₂, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₃, 2*transfinite-1, "Progression",coef)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₃, 4*transfinite-3)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₄, transfinite)


        gmsh.model.mesh.setTransfiniteCurve(Γ₁, 3*transfinite-1)
        gmsh.model.mesh.setTransfiniteCurve(Γ₂, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ₃, 3*transfinite-1)
        gmsh.model.mesh.setTransfiniteCurve(Γ₄, 2*transfinite-1)

        # gmsh.model.mesh.setTransfiniteSurface(Ω)
       
    end

    gmsh.model.mesh.setAlgorithm(2, Ω, 1)
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
    gmsh.model.geo.addPoint(0.0,0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝐿,𝐷, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  𝐿, 𝐷+d, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0, 𝐷, 0.0, lc, 4)
    Γ₁ = gmsh.model.geo.addLine(1, 2, 1)
    Γ₂ = gmsh.model.geo.addLine(2, 3, 2)
    Γ₃ = gmsh.model.geo.addLine(3, 4, 3)
    Γ₄ = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [Γ₁,Γ₃], -1, "Γʳ")
    gmsh.model.addPhysicalGroup(1, [Γ₂], -1, "Γᵗ")
    gmsh.model.addPhysicalGroup(1, [Γ₄], -1, "Γᵍ")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Ω, Γ₁, Γ₂, Γ₃, Γ₄
end

end