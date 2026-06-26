
module Balloon

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 40.0
𝐷 = 10.0
t = 2.0
function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Balloon")

   S1, S2, S3, S4, S5, S6, S7, S8, Γ1, Γ2, Γ3, Γ4,Γ5,Γ6,Γ7,Γ8,Γ9,Γ10,Γ11,Γ12,Γ13,Γ14,Γ15,Γ16, C1, C2, C3, C4,C5,C6,C7,C8= generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    if transfinite > 0
        # gmsh.model.mesh.setTransfiniteCurve(Γ₁, 40*transfinite-3)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₂, 10*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₃, 40*transfinite-3)
        # gmsh.model.mesh.setTransfiniteCurve(Γ₄, 10*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(Γ5, 20*transfinite-3)
        # gmsh.model.mesh.setTransfiniteCurve(Γ6, 6*transfinite)
        # gmsh.model.mesh.setTransfiniteCurve(Γ7, 20*transfinite-3)
        # gmsh.model.mesh.setTransfiniteCurve(Γ8, 6*transfinite)
  
#
        n_long = 20 * transfinite+1
        n_short = 10 * transfinite+1
        n_thick = 2 * transfinite+1
        n_thick2 = 6 * transfinite+1
        # 1. 设置所有线的节点数 (注意：相对的边必须一致)
        for l in [Γ2, Γ13, Γ15, Γ8] 
            gmsh.model.mesh.setTransfiniteCurve(l, n_long) 
        end
        for l in [Γ1, Γ9, Γ3, Γ7, C2,C7,C3,C6] 
            gmsh.model.mesh.setTransfiniteCurve(l, n_short) 
        end
        for l in [C1, C4, C5, C8,Γ12, Γ4, Γ10, Γ6] 
            gmsh.model.mesh.setTransfiniteCurve(l, n_thick) 
        end
        for l in [Γ11, Γ16,Γ14, Γ5,] 
            gmsh.model.mesh.setTransfiniteCurve(l, n_thick2) 
        end




        gmsh.model.mesh.setTransfiniteSurface(S1)
        gmsh.model.mesh.setTransfiniteSurface(S2)
        gmsh.model.mesh.setTransfiniteSurface(S3)
        gmsh.model.mesh.setTransfiniteSurface(S4)
        gmsh.model.mesh.setTransfiniteSurface(S5)
        gmsh.model.mesh.setTransfiniteSurface(S6)
        gmsh.model.mesh.setTransfiniteSurface(S7)
        gmsh.model.mesh.setTransfiniteSurface(S8)
        
       
    end

   
        for S in [S1, S2, S3, S4, S5, S6, S7, S8] 
            gmsh.model.mesh.setAlgorithm(2, S, 1)
        end
    gmsh.model.mesh.generate(2)
    #  gmsh.model.mesh.refine()  
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.setOrder(order)   
    # gmsh.model.mesh.SecondOrderLinear= 1
    # gmsh.model.mesh.secondOrderIncomplete = true
     gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
   gmsh.model.mesh.setOrder(order)


   tag1 = BenchmarkExample.addEdgeElements((2,1), order)
   tag2 = BenchmarkExample.addEdgeElements((2,2), order)
   tag3 = BenchmarkExample.addEdgeElements((2,3), order)
   tag4 = BenchmarkExample.addEdgeElements((2,4), order)
   tag5 = BenchmarkExample.addEdgeElements((2,5), order)
   tag6 = BenchmarkExample.addEdgeElements((2,6), order)
   tag7 = BenchmarkExample.addEdgeElements((2,7), order)
   tag8 = BenchmarkExample.addEdgeElements((2,8), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8], -1, "Γ")
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

    gmsh.model.geo.addPoint(𝐷,-𝐷/2, 0.0, lc, 9)
    gmsh.model.geo.addPoint(  𝐿-𝐷,-𝐷/2, 0.0, lc, 10)
    gmsh.model.geo.addPoint(  𝐿,  -𝐷/2+t, 0.0, lc, 11)
    gmsh.model.geo.addPoint(  𝐿, 𝐷/2-t, 0.0, lc, 12)

    gmsh.model.geo.addPoint(𝐿-𝐷,𝐷/2, 0.0, lc, 13)
    gmsh.model.geo.addPoint(  𝐷,𝐷/2, 0.0, lc, 14)
    gmsh.model.geo.addPoint(  0, 𝐷/2-t, 0.0, lc, 15)
    gmsh.model.geo.addPoint(0.0, -𝐷/2+t, 0.0, lc, 16)

    Γ1 = gmsh.model.geo.addLine(1, 9, 1)
    Γ2 = gmsh.model.geo.addLine(9, 10, 2)
    Γ3 = gmsh.model.geo.addLine(10, 2, 3)
    Γ4 = gmsh.model.geo.addLine(2, 11, 4)
    Γ5 = gmsh.model.geo.addLine(11, 12, 5)
    Γ6 = gmsh.model.geo.addLine(12, 3, 6)
    Γ7 = gmsh.model.geo.addLine(3, 13, 7)
    Γ8 = gmsh.model.geo.addLine(13, 14, 8)
    Γ9 = gmsh.model.geo.addLine(14, 4, 9)
    Γ10 = gmsh.model.geo.addLine(4, 15, 10)
    Γ11 = gmsh.model.geo.addLine(15, 16, 11)
    Γ12 = gmsh.model.geo.addLine(16, 1, 12)

    Γ13 = gmsh.model.geo.addLine(5, 6, 13)
    Γ14 = gmsh.model.geo.addLine(6, 7, 14)
    Γ15 = gmsh.model.geo.addLine(7, 8, 15)
    Γ16 = gmsh.model.geo.addLine(8, 5, 16)



    C1 = gmsh.model.geo.addLine(9, 5, 17)
    C2 = gmsh.model.geo.addLine(5, 16, 18)
    C3 = gmsh.model.geo.addLine(11, 6, 19)
    C4 = gmsh.model.geo.addLine(6, 10, 20)
    C5 = gmsh.model.geo.addLine(13, 7, 21)
    C6 = gmsh.model.geo.addLine(7, 12, 22)
    C7 = gmsh.model.geo.addLine(15, 8, 23)
    C8 = gmsh.model.geo.addLine(8, 14, 24)
   
    # 左下
     gmsh.model.geo.addCurveLoop([Γ12, Γ1, C1, C2],1)
    S1 = gmsh.model.geo.addPlaneSurface([1],1)
    # 下中
    gmsh.model.geo.addCurveLoop([-C1, Γ2, -C4, -Γ13],2)
    S2 = gmsh.model.geo.addPlaneSurface([2],2)
    # 右下
    gmsh.model.geo.addCurveLoop([C3, C4, Γ3, Γ4],3)
    S3 = gmsh.model.geo.addPlaneSurface([3],3)
    # 右中
    gmsh.model.geo.addCurveLoop([-Γ14, -C3, Γ5, -C6],4)
    S4 = gmsh.model.geo.addPlaneSurface([4],4)

    # S5: 右上 (点 12, 3, 13, 7)
     gmsh.model.geo.addCurveLoop([C6, Γ6, Γ7, C5],5)
    S5 = gmsh.model.geo.addPlaneSurface([5],5)

    # S6: 上中 (点 13, 14, 8, 7)
    gmsh.model.geo.addCurveLoop([-C5, Γ8, -C8, -Γ15],6)
    S6 = gmsh.model.geo.addPlaneSurface([6],6)

    # S7: 左上 (点 14, 4, 15, 8)
    gmsh.model.geo.addCurveLoop([C8, Γ9, Γ10, C7],7)
    S7 = gmsh.model.geo.addPlaneSurface([7],7)

    # S8: 左中 (点 15, 16, 5, 8)
    gmsh.model.geo.addCurveLoop([-C7, Γ11, -C2, -Γ16],8)
    S8 = gmsh.model.geo.addPlaneSurface([8],8)


    # gmsh.model.geo.addCurveLoop([1,2,3,4,-5,-6,-7,-8],1)
    
    gmsh.model.geo.synchronize()

    
    gmsh.model.addPhysicalGroup(1, [Γ13], -1, "Γᵗ1")
    gmsh.model.addPhysicalGroup(1, [Γ14], -1, "Γᵗ2")
    gmsh.model.addPhysicalGroup(1, [Γ15], -1, "Γᵗ3")
    gmsh.model.addPhysicalGroup(1, [Γ16], -1, "Γᵗ4")
    gmsh.model.addPhysicalGroup(1, [Γ10], -1, "Γᵍ1")
    gmsh.model.addPhysicalGroup(1, [Γ11], -1, "Γᵍ2")
    gmsh.model.addPhysicalGroup(1, [Γ12], -1, "Γᵍ3")
    gmsh.model.addPhysicalGroup(2, [S1, S2, S3, S4, S5, S6, S7, S8], -1, "Ω")

    return S1, S2, S3, S4, S5, S6, S7, S8, Γ1, Γ2, Γ3, Γ4,Γ5,Γ6,Γ7,Γ8,Γ9,Γ10,Γ11,Γ12,Γ13,Γ14,Γ15,Γ16, C1, C2, C3, C4,C5,C6,C7,C8
end

end