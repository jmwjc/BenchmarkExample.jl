
module Balloon2

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 40.0
𝐷 = 10.0
t = 2.0
function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Balloon2")

  S1,Γ1, Γ2, Γ3, Γ4,Γ5,Γ6,Γ7,Γ8= generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    if transfinite > 0
        gmsh.model.mesh.setTransfiniteCurve(Γ1, 40*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ2, 10*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ3, 40*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ4, 10*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ5, 20*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ6, 6*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ7, 20*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(Γ8, 6*transfinite+1)
  
#

        # gmsh.model.mesh.setTransfiniteSurface(S1)
      
        
       
    end

   
        
    gmsh.model.mesh.setAlgorithm(2, S1, 1)
        
    gmsh.model.mesh.generate(2)
    #  gmsh.model.mesh.refine()  
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.setOrder(order)   
    gmsh.model.mesh.SecondOrderLinear= 1
    # gmsh.model.mesh.secondOrderIncomplete = true
     gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
   gmsh.model.mesh.setOrder(order)


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

    gmsh.model.geo.addPoint(𝐷,-𝐷/2+t, 0.0, lc, 5)
    gmsh.model.geo.addPoint(  𝐿-𝐷,-𝐷/2+t, 0.0, lc, 6)
    gmsh.model.geo.addPoint(  𝐿-𝐷, 𝐷/2-t, 0.0, lc, 7)
    gmsh.model.geo.addPoint(𝐷, 𝐷/2-t, 0.0, lc, 8)

   

    Γ1 = gmsh.model.geo.addLine(1, 2, 1)
    Γ2 = gmsh.model.geo.addLine(2, 3, 2)
    Γ3 = gmsh.model.geo.addLine(3, 4, 3)
    Γ4 = gmsh.model.geo.addLine(4, 1, 4)
    Γ5 = gmsh.model.geo.addLine(5, 6, 5)
    Γ6 = gmsh.model.geo.addLine(6, 7, 6)
    Γ7 = gmsh.model.geo.addLine(7, 8, 7)
    Γ8 = gmsh.model.geo.addLine(8, 5, 8)
   



   
    
     gmsh.model.geo.addCurveLoop([Γ1, Γ2, Γ3, Γ4,-Γ5, -Γ6, -Γ7, -Γ8,],1)
    S1 = gmsh.model.geo.addPlaneSurface([1],1)
   


    # gmsh.model.geo.addCurveLoop([1,2,3,4,-5,-6,-7,-8],1)
    
    gmsh.model.geo.synchronize()

    
    gmsh.model.addPhysicalGroup(1, [Γ5], -1, "Γᵗ1")
    gmsh.model.addPhysicalGroup(1, [Γ6], -1, "Γᵗ2")
    gmsh.model.addPhysicalGroup(1, [Γ7], -1, "Γᵗ3")
    gmsh.model.addPhysicalGroup(1, [Γ8], -1, "Γᵗ4")
    gmsh.model.addPhysicalGroup(1, [Γ4], -1, "Γᵍ")
   
    gmsh.model.addPhysicalGroup(2, [1], -1, "Ω")

    return S1,Γ1, Γ2, Γ3, Γ4,Γ5,Γ6,Γ7,Γ8
end

end