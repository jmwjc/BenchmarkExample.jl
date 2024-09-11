
module PlateWithHole

import ..BenchmarkExample
import Gmsh: gmsh

𝑎 = 1.0
𝑏 = 5.0

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Plate with Hole")

    if transfinite > 0
        Ω₁, Ω₂, Ω₃, Ω₄, Ω₅ = generateGeo(lc, transfinite)
    end

    if quad
        gmsh.model.mesh.setRecombine(2, Ω₁)
        gmsh.model.mesh.setRecombine(2, Ω₂)
        gmsh.model.mesh.setRecombine(2, Ω₃)
        gmsh.model.mesh.setRecombine(2, Ω₄)
        gmsh.model.mesh.setRecombine(2, Ω₅)
    end
    

    gmsh.model.mesh.setAlgorithm(2, Ω₁, 1)
    gmsh.model.mesh.setAlgorithm(2, Ω₂, 1)
    gmsh.model.mesh.setAlgorithm(2, Ω₃, 1)
    gmsh.model.mesh.setAlgorithm(2, Ω₄, 1)
    gmsh.model.mesh.setAlgorithm(2, Ω₅, 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    # gmsh.finalize()
end

@inline function generateGeo(lc = 1.0, n = 1)
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝑎, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  𝑏, 0.0, 0.0, lc, 3)
    gmsh.model.geo.addPoint(  𝑏,   𝑏, 0.0, lc, 4)
    gmsh.model.geo.addPoint(0.0,   𝑏, 0.0, lc, 5)
    gmsh.model.geo.addPoint(0.0,   𝑎, 0.0, lc, 6)
    gmsh.model.geo.addPoint(2*𝑎, 0.0, 0.0, lc, 7)
    gmsh.model.geo.addPoint(0.0, 2*𝑎, 0.0, lc, 8)
    gmsh.model.geo.addPoint(     𝑏, 1.5*𝑎, 0.0, lc, 9)
    gmsh.model.geo.addPoint( 1.5*𝑎,     𝑏, 0.0, lc, 10)
    gmsh.model.geo.addPoint(2^0.5/2*𝑎, 2^0.5/2*𝑎, 0.0, lc, 11)
    gmsh.model.geo.addPoint( 2^0.5*𝑎, 2^0.5*𝑎, 0.0, lc, 12)

    gmsh.model.geo.addLine(2, 7, 1)
    gmsh.model.geo.addLine(7, 3, 2)
    gmsh.model.geo.addLine(3, 9, 3)
    gmsh.model.geo.addLine(9, 4, 4)
    gmsh.model.geo.addLine(4, 10, 5)
    gmsh.model.geo.addLine(10, 5, 6)
    gmsh.model.geo.addLine(5, 8, 7)
    gmsh.model.geo.addLine(8, 6, 8)
    gmsh.model.geo.addCircleArc(6, 1, 11, 9)
    gmsh.model.geo.addCircleArc(11, 1, 2, 10)
    gmsh.model.geo.addCircleArc(7, 1, 12, 11)
    gmsh.model.geo.addCircleArc(12, 1, 8, 12)
    gmsh.model.geo.addLine(11, 12, 13)
    gmsh.model.geo.addLine(12, 9, 14)
    gmsh.model.geo.addLine(10, 12, 15)

    gmsh.model.geo.addCurveLoop([9,13,12,8],1)
    gmsh.model.geo.addCurveLoop([10,1,11,-13],2)
    gmsh.model.geo.addCurveLoop([-12,-15,6,7],3)
    gmsh.model.geo.addCurveLoop([11,14,-3,-2],4)
    gmsh.model.geo.addCurveLoop([14,4,5,15],5)
    Ω₁ = gmsh.model.geo.addPlaneSurface([1],1)
    Ω₂ = gmsh.model.geo.addPlaneSurface([2],2)
    Ω₃ = gmsh.model.geo.addPlaneSurface([3],3)
    Ω₄ = gmsh.model.geo.addPlaneSurface([4],4)
    Ω₅ = gmsh.model.geo.addPlaneSurface([5],5)

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [1,2], -1, "Γᵍ₁")
    gmsh.model.addPhysicalGroup(1, [7,8], -1, "Γᵍ₂")
    gmsh.model.addPhysicalGroup(1, [3,4], -1, "Γᵗ₁")
    gmsh.model.addPhysicalGroup(1, [5,6], -1, "Γᵗ₂")
    gmsh.model.addPhysicalGroup(1, [9,10], -1, "Γᵗ₃")
    gmsh.model.addPhysicalGroup(2, [1,2,3,4,5], -1, "Ω")

    gmsh.model.mesh.setTransfiniteCurve(1, n)
    gmsh.model.mesh.setTransfiniteCurve(3, n)
    gmsh.model.mesh.setTransfiniteCurve(6, n)
    gmsh.model.mesh.setTransfiniteCurve(8, n)
    gmsh.model.mesh.setTransfiniteCurve(9, n)
    gmsh.model.mesh.setTransfiniteCurve(10, n)
    gmsh.model.mesh.setTransfiniteCurve(11, n)
    gmsh.model.mesh.setTransfiniteCurve(12, n)
    gmsh.model.mesh.setTransfiniteCurve(13, n)
    gmsh.model.mesh.setTransfiniteCurve(2, 2*n)
    gmsh.model.mesh.setTransfiniteCurve(4, 2*n)
    gmsh.model.mesh.setTransfiniteCurve(5, 2*n)
    gmsh.model.mesh.setTransfiniteCurve(7, 2*n)
    gmsh.model.mesh.setTransfiniteCurve(14, 2*n)
    gmsh.model.mesh.setTransfiniteCurve(15, 2*n)
    gmsh.model.mesh.setTransfiniteSurface(Ω₁)
    gmsh.model.mesh.setTransfiniteSurface(Ω₂)
    gmsh.model.mesh.setTransfiniteSurface(Ω₃)
    gmsh.model.mesh.setTransfiniteSurface(Ω₄)
    gmsh.model.mesh.setTransfiniteSurface(Ω₅)

    return Ω₁, Ω₂, Ω₃, Ω₄, Ω₅
end
end