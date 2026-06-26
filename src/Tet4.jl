
module Tet4

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 1.0


function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
   
    gmsh.model.add("Tet4")
    Ω,  Γ1, Γ2, Γ3, Γ4,L1,L2,L3,L4,L5,L6= generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(3, Ω)
    end
    
    if transfinite > 0
        # gmsh.model.mesh.setTransfiniteCurve(L1, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L2, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L3, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L4, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L5, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L6, transfinite)
      
       


       
        # gmsh.model.mesh.setTransfiniteSurface(Γ1)
        # gmsh.model.mesh.setTransfiniteSurface(Γ2)
        # gmsh.model.mesh.setTransfiniteSurface(Γ3)
        # gmsh.model.mesh.setTransfiniteSurface(Γ4)
        
        
        # gmsh.model.mesh.setTransfiniteVolume(Ω)
       
        

    end
    
    
    # gmsh.model.geo.synchronize()
    gmsh.model.mesh.setAlgorithm(3, Ω, 1)
    # gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    # gmsh.model.mesh.SecondOrderLinear = 1
    # gmsh.model.mesh.secondOrderIncomplete = true
    # gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.optimize("Netgen")
    tag = BenchmarkExample.addFaceElements((3,1), order)
    gmsh.model.addPhysicalGroup(2, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
   

    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    gmsh.model.geo.addPoint(0.0,0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝐿,0.0, 0.0, lc, 2)
   
    gmsh.model.geo.addPoint(0.0,  𝐿, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0,  0.0, 𝐿, lc, 4)
   


    L1  = gmsh.model.geo.addLine( 1,  2, 1) 
    L2  = gmsh.model.geo.addLine( 2,  3, 2) 
    L3  = gmsh.model.geo.addLine( 3,  1, 3) 
    L4  = gmsh.model.geo.addLine( 1, 4, 4) 
    L5  = gmsh.model.geo.addLine( 2, 4, 5) 
    L6  = gmsh.model.geo.addLine( 3, 4, 6) 
  

    gmsh.model.geo.addCurveLoop([1,2,3],1)
    gmsh.model.geo.addCurveLoop([1,5,-4],2)
    gmsh.model.geo.addCurveLoop([2,6,-5],3)
    gmsh.model.geo.addCurveLoop([3,4,-6,],4)
    
    

    Γ1 = gmsh.model.geo.addPlaneSurface([1],1)
    Γ2 = gmsh.model.geo.addPlaneSurface([2],2)
    Γ3 = gmsh.model.geo.addPlaneSurface([3],3)
    Γ4 = gmsh.model.geo.addPlaneSurface([4],4)
   
   

    gmsh.model.geo.addSurfaceLoop([1,2,3,4],1)
    
    Ω = gmsh.model.geo.addVolume([1],1)
   
    gmsh.model.geo.synchronize()

    
    gmsh.model.addPhysicalGroup(2, [Γ1,Γ2,Γ3,Γ4], -1, "Γᵍ")
    
    
    gmsh.model.addPhysicalGroup(3, [Ω], -1, "Ω")

    return Ω,  Γ1, Γ2, Γ3, Γ4,L1,L2,L3,L4,L5,L6
end
end