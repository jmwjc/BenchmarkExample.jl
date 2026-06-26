
module TickWalledCylinder

import ..BenchmarkExample
import Gmsh: gmsh
# import Statistics
r_in= 3.0
r_out= 9.0
deg=90*π /180 
function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false, mode = 1, coef = 1.0)
    gmsh.initialize()
    gmsh.model.add("TickWalledCylinder")
    Ω, Γ1, Γ2, Γ3, Γ4 = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    if transfinite > 0
      
        gmsh.model.mesh.setTransfiniteCurve(Γ1, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ2, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ3, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(Γ4, 3*transfinite)

        # gmsh.model.mesh.setTransfiniteCurve(Γ1, 6*transfinite, "Progression",coef)
        # gmsh.model.mesh.setTransfiniteCurve(Γ2, 6*transfinite, "Progression",-coef)
        # gmsh.model.mesh.setTransfiniteCurve(Γ3, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(Γ4, transfinite)
        # gmsh.model.mesh.setTransfiniteSurface(Ω)

        # gmsh.model.mesh.setTransfiniteCurve(Γ1, transfinite, "Progression",coef)
        # gmsh.model.mesh.setTransfiniteCurve(Γ2, transfinite, "Progression",-coef)
        # gmsh.model.mesh.setTransfiniteCurve(Γ3, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(Γ4, transfinite)
        # gmsh.model.mesh.setTransfiniteSurface(Ω)
       
       
    end
    gmsh.model.mesh.setAlgorithm(2, Ω, 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    gmsh.model.mesh.SecondOrderLinear = 1
    gmsh.model.mesh.secondOrderIncomplete = true
    #  gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
   
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")

    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc=1.0)
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint( r_out, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(r_out * cos(deg), r_out * sin(deg), 0, lc, 3)
    gmsh.model.geo.addPoint( r_in, 0, 0, lc, 4)
    gmsh.model.geo.addPoint(r_in * cos(deg), r_in * sin(deg), 0, lc, 5)
   
    Γ1=gmsh.model.geo.addLine(4, 2, 1)
    Γ2=gmsh.model.geo.addLine(3, 5, 2)
    Γ3=gmsh.model.geo.addCircleArc(5, 1, 4, 3)
    Γ4=gmsh.model.geo.addCircleArc(2, 1, 3, 4)
   

    gmsh.model.geo.addCurveLoop([1,4,2,3],1)
    
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
  

    gmsh.model.geo.synchronize()

    
        gmsh.model.addPhysicalGroup(1, [Γ1], -1, "Γ1")
        gmsh.model.addPhysicalGroup(1, [Γ2], -1, "Γ2")
        gmsh.model.addPhysicalGroup(1, [Γ3], -1, "Γ3")
        gmsh.model.addPhysicalGroup(1, [Γ4], -1, "Γ4")
        gmsh.model.addPhysicalGroup(2, [1], -1, "Ω")
    


    return  Ω, Γ1, Γ2, Γ3, Γ4
end

end