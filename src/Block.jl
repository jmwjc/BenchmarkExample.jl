
module Block

import ..BenchmarkExample
import Gmsh: gmsh

𝐿 = 1.0


function generateMsh(filepath::String; lc = 1.0, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Block")

    Ω1, Ω2, Ω3, Ω4, Γ1, Γ2, Γ3, Γ4, Γ5, Γ6, Γ7, Γ8, Γ9, Γ10, Γ11, Γ12, Γ13, Γ14, Γ15, Γ16, Γ17, Γ18, Γ19, Γ20, L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,L21,L22,L23,L24,L25,L26,L27,L28,L29,L30,L31,L32,L33 = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(3, Ω1)
        gmsh.model.mesh.setRecombine(3, Ω2)
        gmsh.model.mesh.setRecombine(3, Ω3)
        gmsh.model.mesh.setRecombine(3, Ω4)
        # gmsh.model.mesh.setRecombine(3, Ω)
    end
    
    if transfinite > 0
        # gmsh.model.mesh.setTransfiniteCurve([L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,L21,L22,L23,L24], transfinite)
        gmsh.model.mesh.setTransfiniteCurve(L1, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L2, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L3, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L4, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L5, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L6, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L7, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L8, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L9, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L10, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L11, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L12, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L13, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L14, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L15, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L16, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L17, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L18, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L19, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L20, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L21, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L22, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L23, transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L24, transfinite+1)


        gmsh.model.mesh.setTransfiniteCurve(L25, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L26, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L27, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L28, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L29, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L30, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L31, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L32, 2*transfinite+1)
        gmsh.model.mesh.setTransfiniteCurve(L33, 2*transfinite+1)
        
        
        
        
        
        # gmsh.model.mesh.setTransfiniteSurface([Γ1, Γ2, Γ3, Γ4, Γ5, Γ6, Γ7, Γ8, Γ9, Γ10, Γ11, Γ12, Γ13, Γ14, Γ15, Γ16, Γ17, Γ18, Γ19, Γ20])
        gmsh.model.mesh.setTransfiniteSurface(Γ1)
        gmsh.model.mesh.setTransfiniteSurface(Γ1)
        gmsh.model.mesh.setTransfiniteSurface(Γ2)
        gmsh.model.mesh.setTransfiniteSurface(Γ3)
        gmsh.model.mesh.setTransfiniteSurface(Γ4)
        # gmsh.model.mesh.setTransfiniteSurface(Γ5)
        # gmsh.model.mesh.setTransfiniteSurface(Γ6)
        # gmsh.model.mesh.setTransfiniteSurface(Γ7)
        # gmsh.model.mesh.setTransfiniteSurface(Γ8)
        gmsh.model.mesh.setTransfiniteSurface(Γ5, "Right",  [5, 14, 18, 17])
        gmsh.model.mesh.setTransfiniteSurface(Γ6, "Left", [14, 6, 15, 18])
        gmsh.model.mesh.setTransfiniteSurface(Γ7, "Right",  [15, 7, 16, 18])
        gmsh.model.mesh.setTransfiniteSurface(Γ8, "Left", [16, 8, 17, 18])
        gmsh.model.mesh.setTransfiniteSurface(Γ9)
        gmsh.model.mesh.setTransfiniteSurface(Γ10)
        gmsh.model.mesh.setTransfiniteSurface(Γ11)
        gmsh.model.mesh.setTransfiniteSurface(Γ12)
        # gmsh.model.mesh.setTransfiniteSurface(Γ13)
        gmsh.model.mesh.setTransfiniteSurface(Γ13, "Right", [3, 11, 16, 7])
        gmsh.model.mesh.setTransfiniteSurface(Γ14)
        gmsh.model.mesh.setTransfiniteSurface(Γ15)
        gmsh.model.mesh.setTransfiniteSurface(Γ16)
        gmsh.model.mesh.setTransfiniteSurface(Γ17)
        gmsh.model.mesh.setTransfiniteSurface(Γ18)
        gmsh.model.mesh.setTransfiniteSurface(Γ19)
        gmsh.model.mesh.setTransfiniteSurface(Γ20)

        gmsh.model.mesh.setTransfiniteVolume(Ω1)
        gmsh.model.mesh.setTransfiniteVolume(Ω2)
        gmsh.model.mesh.setTransfiniteVolume(Ω3)
        gmsh.model.mesh.setTransfiniteVolume(Ω4)
        

    end

    gmsh.model.mesh.setAlgorithm(3, Ω1, 1)
    gmsh.model.mesh.setAlgorithm(3, Ω2, 1)
    gmsh.model.mesh.setAlgorithm(3, Ω3, 1)
    gmsh.model.mesh.setAlgorithm(3, Ω4, 1)
   
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    # gmsh.model.mesh.SecondOrderLinear = 1
    # gmsh.model.mesh.secondOrderIncomplete = true

    
    gmsh.option.setNumber("Mesh.Optimize", 1)
   
    tag1 = BenchmarkExample.addFaceElements((3,1), order)
    tag2 = BenchmarkExample.addFaceElements((3,2), order)
    tag3 = BenchmarkExample.addFaceElements((3,3), order)
    tag4 = BenchmarkExample.addFaceElements((3,4), order)
    gmsh.model.addPhysicalGroup(2, [tag1,tag2,tag3,tag4], -1, "Γ")
    gmsh.model.geo.synchronize()
    
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    gmsh.model.geo.addPoint(0.0,0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  𝐿,0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  𝐿,  𝐿, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0,  𝐿, 0.0, lc, 4)
    gmsh.model.geo.addPoint(0.0,  0.0, 𝐿, lc, 5)
    gmsh.model.geo.addPoint(  𝐿,  0.0, 𝐿, lc, 6)
    gmsh.model.geo.addPoint(  𝐿,    𝐿, 𝐿, lc, 7)
    gmsh.model.geo.addPoint(0.0,    𝐿, 𝐿, lc, 8)
    gmsh.model.geo.addPoint(0.5*𝐿,   0.0,   0.0, lc, 9)
    gmsh.model.geo.addPoint(    𝐿, 0.5*𝐿,   0.0, lc, 10)
    gmsh.model.geo.addPoint(0.5*𝐿,     𝐿,   0.0, lc, 11)
    gmsh.model.geo.addPoint(  0.0, 0.5*𝐿,   0.0, lc, 12)
    gmsh.model.geo.addPoint(0.5*𝐿, 0.5*𝐿,   0.0, lc, 13)
    gmsh.model.geo.addPoint(0.5*𝐿,   0.0,     𝐿, lc, 14)
    gmsh.model.geo.addPoint(    𝐿, 0.5*𝐿,     𝐿, lc, 15)
    gmsh.model.geo.addPoint(0.5*𝐿,     𝐿,     𝐿, lc, 16)
    gmsh.model.geo.addPoint(  0.0, 0.5*𝐿,     𝐿, lc, 17)
    gmsh.model.geo.addPoint(0.5*𝐿, 0.5*𝐿,     𝐿, lc, 18)


    L1  = gmsh.model.geo.addLine( 1,  9, 1) 
    L2  = gmsh.model.geo.addLine( 9,  2, 2) 
    L3  = gmsh.model.geo.addLine( 2, 10, 3) 
    L4  = gmsh.model.geo.addLine( 10, 3, 4) 
    L5  = gmsh.model.geo.addLine( 3, 11, 5) 
    L6  = gmsh.model.geo.addLine(11, 4 , 6) 
    L7  = gmsh.model.geo.addLine( 4, 12, 7) 
    L8  = gmsh.model.geo.addLine(12,  1, 8) 
    L9  = gmsh.model.geo.addLine( 9,  13, 9)
    L10 = gmsh.model.geo.addLine(13, 12, 10)  
    L11 = gmsh.model.geo.addLine(11, 13, 11) 
    L12 = gmsh.model.geo.addLine(13, 10, 12)
    L13 = gmsh.model.geo.addLine( 5, 14, 13)
    L14 = gmsh.model.geo.addLine(14,  6, 14)
    L15 = gmsh.model.geo.addLine( 6, 15, 15)
    L16 = gmsh.model.geo.addLine(15,  7, 16)
    L17 = gmsh.model.geo.addLine( 7, 16, 17)
    L18 = gmsh.model.geo.addLine(16,  8, 18)
    L19 = gmsh.model.geo.addLine( 8, 17, 19)
    L20 = gmsh.model.geo.addLine(17,  5, 20)
    L21 = gmsh.model.geo.addLine(14, 18, 21)
    L22 = gmsh.model.geo.addLine(18, 17, 22)
    L23 = gmsh.model.geo.addLine(16, 18, 23)
    L24 = gmsh.model.geo.addLine(18, 15, 24)
    L25 = gmsh.model.geo.addLine( 1,  5, 25)
    L26 = gmsh.model.geo.addLine( 2,  6, 26)
    L27 = gmsh.model.geo.addLine( 3,  7, 27)
    L28 = gmsh.model.geo.addLine( 4,  8, 28)
    L29 = gmsh.model.geo.addLine( 9, 14, 29)
    L30 = gmsh.model.geo.addLine(10, 15, 30)
    L31 = gmsh.model.geo.addLine(11, 16, 31)
    L32 = gmsh.model.geo.addLine(12, 17, 32)
    L33 = gmsh.model.geo.addLine(13, 18, 33)

    gmsh.model.geo.addCurveLoop([1,9,10,8],1)
    gmsh.model.geo.addCurveLoop([2,3,-12,-9],2)
    gmsh.model.geo.addCurveLoop([4,5,11,12],3)
    gmsh.model.geo.addCurveLoop([6,7,-10,-11],4)
    gmsh.model.geo.addCurveLoop([13,21,22,20],5)
    gmsh.model.geo.addCurveLoop([14,15,-24,-21],6)
    gmsh.model.geo.addCurveLoop([16,17,23,24],7)
    gmsh.model.geo.addCurveLoop([18,19,-22,-23],8)
    gmsh.model.geo.addCurveLoop([1,29,-13,-25],9)
    gmsh.model.geo.addCurveLoop([2,26,-14,-29],10)
    gmsh.model.geo.addCurveLoop([3,30,-15,-26],11)
    gmsh.model.geo.addCurveLoop([4,27,-16,-30],12)
    gmsh.model.geo.addCurveLoop([5,31,-17,-27],13)
    gmsh.model.geo.addCurveLoop([6,28,-18,-31],14)
    gmsh.model.geo.addCurveLoop([7,32,-19,-28],15)
    gmsh.model.geo.addCurveLoop([8,25,-20,-32],16)
    gmsh.model.geo.addCurveLoop([9,33,-21,-29],17)
    gmsh.model.geo.addCurveLoop([10,32,-22,-33],18)
    gmsh.model.geo.addCurveLoop([11,33,-23,-31],19)
    gmsh.model.geo.addCurveLoop([12,30,-24,-33],20)

    Γ1 = gmsh.model.geo.addPlaneSurface([-1],1)
    Γ2 = gmsh.model.geo.addPlaneSurface([-2],2)
    Γ3 = gmsh.model.geo.addPlaneSurface([-3], 3)
    Γ4 = gmsh.model.geo.addPlaneSurface([-4], 4)
    Γ5 = gmsh.model.geo.addPlaneSurface([5], 5)
    Γ6 = gmsh.model.geo.addPlaneSurface([6], 6)
    Γ7 = gmsh.model.geo.addPlaneSurface([7], 7)
    Γ8 = gmsh.model.geo.addPlaneSurface([8], 8)
    Γ9 = gmsh.model.geo.addPlaneSurface([9], 9)
    Γ10 = gmsh.model.geo.addPlaneSurface([10], 10)
    Γ11 = gmsh.model.geo.addPlaneSurface([11], 11)
    Γ12 = gmsh.model.geo.addPlaneSurface([12], 12)
    Γ13 = gmsh.model.geo.addPlaneSurface([13], 13)
    Γ14 = gmsh.model.geo.addPlaneSurface([14], 14)
    Γ15 = gmsh.model.geo.addPlaneSurface([15], 15)
    Γ16 = gmsh.model.geo.addPlaneSurface([16], 16)
    Γ17 = gmsh.model.geo.addPlaneSurface([17], 17)
    Γ18 = gmsh.model.geo.addPlaneSurface([18], 18)
    Γ19 = gmsh.model.geo.addPlaneSurface([19], 19)
    Γ20 = gmsh.model.geo.addPlaneSurface([20], 20)

    gmsh.model.geo.addSurfaceLoop([1,5,9,16,17,18],1 )                                                                                                                                                                                                     
    gmsh.model.geo.addSurfaceLoop([2,6,10,11,-20,-17],2)
    gmsh.model.geo.addSurfaceLoop([3,7,12,13,19,20],3)
    gmsh.model.geo.addSurfaceLoop([4,8,14,15,-18,-19],4)

    Ω1 = gmsh.model.geo.addVolume([1],1)
    Ω2 = gmsh.model.geo.addVolume([2],2)
    Ω3 = gmsh.model.geo.addVolume([3],3)
    Ω4 = gmsh.model.geo.addVolume([4],4)

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [Γ5], -1, "Γᵗ")
    gmsh.model.addPhysicalGroup(2, [Γ1,Γ2,Γ3,Γ4,Γ6,Γ7,Γ8,Γ9,Γ10,Γ15,Γ16], -1, "Γᵍ")
    gmsh.model.addPhysicalGroup(2, [Γ11,Γ12,Γ13,Γ14], -1, "Γʳ")
    

    gmsh.model.addPhysicalGroup(3, [Ω1,Ω2,Ω3,Ω4], -1, "Ω")

    return Ω1, Ω2, Ω3, Ω4, Γ1, Γ2, Γ3, Γ4, Γ5, Γ6, Γ7, Γ8, Γ9, Γ10, Γ11, Γ12, Γ13, Γ14, Γ15, Γ16, Γ17, Γ18, Γ19, Γ20, L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,L21,L22,L23,L24,L25,L26,L27,L28,L29,L30,L31,L32,L33
end
end