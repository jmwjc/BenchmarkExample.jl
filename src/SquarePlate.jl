module SquarePlate
    
import ..BenchmarkExample
import Gmsh: gmsh

const 𝐸 = 10.92e6
const 𝜈 = 0.3
const 𝐿 = 1.0
const ℎ = 0.1
# const 𝑤(x,y) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
# const 𝜃₁(x,y) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
# const 𝜃₂(x,y) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
# const 𝐹(x,y) = E/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Square plate")

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
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    # gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    𝐿 = 1.0
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(𝐿, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(𝐿,  𝐿, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0, 𝐿, 0.0, lc, 4)
    Γᵇ = gmsh.model.geo.addLine(1, 2, 1)
    Γʳ = gmsh.model.geo.addLine(2, 3, 2)
    Γᵗ = gmsh.model.geo.addLine(3, 4, 3)
    Γˡ = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [Γᵇ], -1, "Γᵇ")
    gmsh.model.addPhysicalGroup(1, [Γᵗ], -1, "Γᵗ")
    gmsh.model.addPhysicalGroup(1, [Γˡ], -1, "Γˡ")
    gmsh.model.addPhysicalGroup(1, [Γʳ], -1, "Γʳ")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Γᵇ, Γʳ, Γᵗ, Γˡ, Ω
end

end