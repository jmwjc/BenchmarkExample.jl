
module Isolation_half

import ..BenchmarkExample
import Gmsh: gmsh
W = 2.0
𝐿 = 1.0
H = 1.0

function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
   
    gmsh.model.add("Isolation_half")
    Ω,  Γ1, Γ2, Γ3, Γ4, Γ5, Γ6,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12= generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(3, Ω)
        gmsh.model.mesh.setRecombine(2, Γ1)
        gmsh.model.mesh.setRecombine(2, Γ2)
        gmsh.model.mesh.setRecombine(2, Γ3)
        gmsh.model.mesh.setRecombine(2, Γ4)
        gmsh.model.mesh.setRecombine(2, Γ5)
        gmsh.model.mesh.setRecombine(2, Γ6)
    end
    
    if transfinite > 0
        # gmsh.model.mesh.setTransfiniteCurve(L1, 2*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L2, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L3, 2*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L4, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L5, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L6, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L7, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L8, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L9, 2*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L10, transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L11, 2*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(L12, transfinite)
       

       
        gmsh.model.mesh.setTransfiniteCurve(L1, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L2, 2*transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L3, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L4, 2*transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L5, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L6, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L7, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L8, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L9, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L10, 2*transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L11, transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L12, 2*transfinite)


        gmsh.model.mesh.setTransfiniteSurface(Γ1)
        gmsh.model.mesh.setTransfiniteSurface(Γ2)
        gmsh.model.mesh.setTransfiniteSurface(Γ3)
        gmsh.model.mesh.setTransfiniteSurface(Γ4)
        gmsh.model.mesh.setTransfiniteSurface(Γ5)
        gmsh.model.mesh.setTransfiniteSurface(Γ6)
        
        gmsh.model.mesh.setTransfiniteVolume(Ω)
       
        

    end
    
    
    # gmsh.model.geo.synchronize()
    gmsh.model.mesh.setAlgorithm(3, Ω, 1)
    # gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    # gmsh.model.mesh.SecondOrderLinear = 1
    # gmsh.model.mesh.secondOrderIncomplete = true
    gmsh.option.setNumber("Mesh.Optimize", 1)
    # gmsh.model.mesh.optimize("Netgen")

    tag = BenchmarkExample.addFaceElements((3,1), order)
    gmsh.model.addPhysicalGroup(2, [tag], -1, "Γ")

    gmsh.model.geo.synchronize()
   

    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    gmsh.model.geo.addPoint(0.0,W, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝐿,W, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  𝐿,  2*W, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0,  2*W, 0.0, lc, 4)
    gmsh.model.geo.addPoint(0.0,  W, H, lc, 5)
    gmsh.model.geo.addPoint(  𝐿,  W, H, lc, 6)
    gmsh.model.geo.addPoint(  𝐿,    2*W, H, lc, 7)
    gmsh.model.geo.addPoint(0.0,    2*W, H, lc, 8)
   


    L1  = gmsh.model.geo.addLine( 1,  2, 1) 
    L2  = gmsh.model.geo.addLine( 2,  3, 2) 
    L3  = gmsh.model.geo.addLine( 3, 4, 3) 
    L4  = gmsh.model.geo.addLine( 4, 1, 4) 
    L5  = gmsh.model.geo.addLine( 1, 5, 5) 
    L6  = gmsh.model.geo.addLine( 2, 6, 6) 
    L7  = gmsh.model.geo.addLine( 3, 7, 7) 
    L8  = gmsh.model.geo.addLine(4,  8, 8) 
    L9  = gmsh.model.geo.addLine( 5,  6, 9)
    L10 = gmsh.model.geo.addLine(6, 7, 10)  
    L11 = gmsh.model.geo.addLine(7, 8, 11) 
    L12 = gmsh.model.geo.addLine(8, 5, 12)
    

     gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1),   # 底面 (1→2→3→4)
    gmsh.model.geo.addCurveLoop([5, 9, -6, -1], 2), # 前面 (5→9→-6→-1)
    gmsh.model.geo.addCurveLoop([6, 10, -7, -2], 3),# 右面 (6→10→-7→-2)
    gmsh.model.geo.addCurveLoop([7, 11, -8, -3], 4),# 后面 (7→11→-8→-3)
    gmsh.model.geo.addCurveLoop([8, 12, -5, -4], 5),# 左面 (8→12→-5→-4)
    gmsh.model.geo.addCurveLoop([9, 10, 11, 12], 6) # 顶面 (9→10→11→12)
    

    # gmsh.model.geo.addCurveLoop([-4, -3, -2, -1], 1) # 底面：改为顺时针，法向向下(朝外)
    # gmsh.model.geo.addCurveLoop([1, 6, -9, -5],   2) # 前面：调整为朝前(朝外)
    # gmsh.model.geo.addCurveLoop([2, 7, -10, -6],  3) # 右面：调整为朝右(朝外)
    # gmsh.model.geo.addCurveLoop([3, 8, -11, -7],  4) # 后面：调整为朝后(朝外)
    # gmsh.model.geo.addCurveLoop([4, 5, -12, -8],  5) # 左面：调整为朝左(朝外)
    # gmsh.model.geo.addCurveLoop([9, 10, 11, 12],  6) # 顶面：保持原样，法向向上(朝外)
    Γ1 = gmsh.model.geo.addPlaneSurface([-1],1)
    Γ2 = gmsh.model.geo.addPlaneSurface([-2],2)
    Γ3 = gmsh.model.geo.addPlaneSurface([-3],3)
    Γ4 = gmsh.model.geo.addPlaneSurface([-4],4)
    Γ5 = gmsh.model.geo.addPlaneSurface([-5],5)
    Γ6 = gmsh.model.geo.addPlaneSurface([6],6)
   

    gmsh.model.geo.addSurfaceLoop([1,2,3,4,5,6],1)
    
    Ω = gmsh.model.geo.addVolume([1],1)
   
    gmsh.model.geo.synchronize()

    
    # gmsh.model.addPhysicalGroup(2, [Γ1,Γ2,Γ3,Γ4,Γ5,Γ6], -1, "Γᵍ")
    gmsh.model.addPhysicalGroup(2, [Γ1], -1, "Γᵍ₁")
    gmsh.model.addPhysicalGroup(2, [Γ2], -1, "Γ2")
    gmsh.model.addPhysicalGroup(2, [Γ3], -1, "Γ3")
    gmsh.model.addPhysicalGroup(2, [Γ4], -1, "Γ4")
    gmsh.model.addPhysicalGroup(2, [Γ5], -1, "Γ5")
    gmsh.model.addPhysicalGroup(2, [Γ6], -1, "Γᵍ₂")
    gmsh.model.addPhysicalGroup(3, [Ω], -1, "Ω")

    return Ω,  Γ1, Γ2, Γ3, Γ4, Γ5, Γ6,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12
end
end