module Circular
    
import ..BenchmarkExample
import Gmsh: gmsh

const 𝐸 = 10.92
const 𝜈 = 0.3
const 𝑅 = 5.0
const ℎ = 0.1  # 0.1/1.0/2.5
const 𝑤 = 0.0
const 𝜃₁ = 0.0
const 𝜃₂ = 0.0
const 𝐹 = 1.0
                   #0.1         1.0           2.5       
const 𝑣ᴱ = 39831.0   #𝑣ᴰ = 39831  𝑣ᴰ = 41.599   𝑣ᴰ = 3.262

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Circular")

    generateGeo(lc)

    if transfinite > 0
        for i in 1:9
            gmsh.model.mesh.setTransfiniteCurve(i, transfinite)
        end
        for i in 1:3
            i == 2 ? gmsh.model.mesh.setTransfiniteSurface(i,"Right") : gmsh.model.mesh.setTransfiniteSurface(i)
        end
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
    𝑅 = 5.0
    𝐴 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(𝑅/2, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(𝑅, 0.0, 0.0, lc, 3)
    gmsh.model.geo.addPoint(𝑅/2^0.5, 𝑅/2^0.5, 0.0, lc, 4)
    gmsh.model.geo.addPoint(0.0, 𝑅, 0.0, lc, 5)
    gmsh.model.geo.addPoint(0.0, 𝑅/2, 0.0, lc, 6)
    gmsh.model.geo.addPoint(2.0, 2.0, 0.0, lc, 7)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addCircleArc(3, 1, 4, 3)
    gmsh.model.geo.addCircleArc(4, 1, 5, 4)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(6, 1, 6)
    gmsh.model.geo.addLine(7, 2, 7)
    gmsh.model.geo.addLine(7, 6, 8)
    gmsh.model.geo.addLine(7, 4, 9)
    gmsh.model.geo.addCurveLoop([1,-7,8,6],1)
    gmsh.model.geo.addCurveLoop([2,3,-9,7],2)
    gmsh.model.geo.addCurveLoop([9,4,5,-8],3)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.addPlaneSurface([2],2)
    gmsh.model.geo.addPlaneSurface([3],3)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(0, [𝐴], -1, "𝐴")
    gmsh.model.addPhysicalGroup(1, [1,2], -1, "Γᵇ")
    gmsh.model.addPhysicalGroup(1, [3,4], -1, "Γᵉ")
    gmsh.model.addPhysicalGroup(1, [5,6], -1, "Γˡ")
    gmsh.model.addPhysicalGroup(2, [1,2,3], -1, "Ω")
end

end