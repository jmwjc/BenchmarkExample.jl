module PatchTestThinShell

import ..BenchmarkExample
import Gmsh: gmsh

𝐸 = 1.0
𝜈 = 0.0
𝑅 = 1.0
𝐿 = 0.5*π
𝜃 = 0.5*π
ℎ = 1.0

function variables(cs::Function,𝒖::Function)
    ∇𝒖(x) = gradient(𝒖,x)
    ∇²𝒖(x) = gradient(∇𝒖,x)
    function 𝜺(x)
        grad𝒖 = ∇𝒖(x)
        ∂₁𝒖 = Tensor{1,3}((grad𝒖[:,1]))
        ∂₂𝒖 = Tensor{1,3}((grad𝒖[:,2]))
        ε₁₁ = cs.𝒂₍₁₎(x)⋅∂₁𝒖
        ε₂₂ = cs.𝒂₍₂₎(x)⋅∂₂𝒖
        ε₁₂ = 0.5*(cs.𝒂₍₁₎(x)⋅∂₂𝒖 + cs.𝒂₍₂₎(x)⋅∂₁𝒖)
        return Tensor{2,2}((ε₁₁,ε₁₂,ε₁₂,ε₂₂))
    end
    function 𝜿(x)
        grad𝒖 = ∇𝒖(x)
        grad2𝒖 = ∇²𝒖(x)
        ∂₁𝒖 = Tensor{1,3}((grad𝒖[:,1]))
        ∂₂𝒖 = Tensor{1,3}((grad𝒖[:,2]))
        ∂₁₁𝒖 = Tensor{1,3}((grad2𝒖[:,1,1]))
        ∂₂₂𝒖 = Tensor{1,3}((grad2𝒖[:,2,2]))
        ∂₁₂𝒖 = Tensor{1,3}((grad2𝒖[:,1,2]))
        Γ¹₁₁ = cs.Γ¹₁₁(x)
        Γ¹₁₂ = cs.Γ¹₁₂(x)
        Γ²₂₁ = cs.Γ²₁₂(x)
        Γ²₂₂ = cs.Γ²₂₂(x)
        κ₁₁ = (Γ¹₁₁.*∂₁𝒖+Γ²₁₁.*∂₂𝒖 - ∂₁₁𝒖)⋅𝒂₃(x)
        κ₂₂ = (Γ¹₂₂.*∂₁𝒖+Γ²₂₂.*∂₂𝒖 - ∂₂₂𝒖)⋅𝒂₃(x)
        κ₁₂ = (Γ¹₁₂.*∂₁𝒖+Γ²₁₂.*∂₂𝒖 - ∂₁₂𝒖)⋅𝒂₃(x)
        return Tensor{2,2}((κ₁₁,κ₁₂,κ₁₂,κ₂₂))
    end
    function 𝑪(x)
        a¹¹ = cs.a¹¹(x)
        a²² = cs.a²²(x)
        a¹² = cs.a¹²(x)
        C¹¹¹¹ = a¹¹*a¹¹*Dᵇ
        C²²²² = a²²*a²²*Dᵇ
        C¹¹²² = (ν*a¹¹*a²² + (1-ν)*a¹²*a¹²)*Dᵇ
        C¹¹¹² = a¹¹*a¹²*Dᵇ
        C²²¹² = a²²*a¹²*Dᵇ
        C¹²¹² = (0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²))*Dᵇ
        return SymmetricTensor{4,2}((C¹¹¹¹,C¹¹¹²,C¹¹²²,C¹¹¹²,C¹²¹²,C²²¹²,C¹¹²²,C²²¹²,C²²²²))
    end
    function 𝑵(x)
        Dᵐ = 𝐸*ℎ
        ε = 𝜺(x)
        return Tensor{2,2}((Dᵐ.*(𝑪 ⊡ ε)))
    end
    function 𝑴(x)
        Dᵇ = 𝐸*ℎ^3/12
        κ = 𝜿(x)
        return Tensor{2,2}((Dᵇ.*(𝑪 ⊡ κ)))
    end
    function 𝜃ₙ(x)
        n¹ = cs.n¹(x)
        n² = cs.n²(x)
        𝒂₃ = cs.𝒂₃(x)
        grad𝒖 = ∇𝒖(x)
        ∂₁𝒖 = Tensor{1,3}((grad𝒖[:,1]))
        ∂₂𝒖 = Tensor{1,3}((grad𝒖[:,2]))
        return (n¹.*∂₁𝒖 + n².*∂₂𝒖)⋅𝒂₃
    end
    function 𝒃(x)
        𝒂₁ = cs.𝒂₁(x)
        𝒂₂ = cs.𝒂₂(x)
        𝒂₃ = cs.𝒂₃(x)
        b₁₁ = cs.b₁₁(x)
        b₂₂ = cs.b₂₂(x)
        b₁₂ = cs.b₁₂(x)
        Γ¹₁₁ = cs.Γ¹₁₁(x)
        Γ¹₁₂ = cs.Γ¹₁₂(x)
        Γ¹₂₂ = cs.Γ¹₂₂(x)
        Γ²₁₁ = cs.Γ²₁₁(x)
        Γ²₂₁ = cs.Γ²₁₂(x)
        Γ²₂₂ = cs.Γ²₂₂(x)
        ∂₁Γ¹₁₁ = cs.∂₁Γ¹₁₁(x)
        ∂₂Γ¹₁₁ = cs.∂₂Γ¹₁₁(x)
        Γ¹ = Tensor{2,2}((Γ¹₁₁,Γ¹₁₂,Γ¹₁₂,Γ¹₂₂))
        Γ² = Tensor{2,2}((Γ²₁₁,Γ²₂₁,Γ²₂₁,Γ²₂₂))
        b = Tensor{2,2}((b₁₁,b₁₂,b₁₂,b₂₂))
        Γᵞᵧ = Vec{2}((Γ¹₁₁+Γ²₂₁,Γ¹₁₂+Γ²₂₂))
        N = 𝑵(x)
        ∇N = gradient(𝑵,x)
        divN = ∇N[:,1,1] + ∇N[:,2,2]
        NΓᵞᵧ = N ⋅ Γᵞᵧ
        𝒃1 = (Γ¹ ⊡ N + divN[1] + NΓᵞᵧ[1]).*𝒂₁ + (Γ² ⊡ N + divN[2] + NΓᵞᵧ[2]).*𝒂₂ + b ⊡ N.*𝒂₃
    end
    return ()->(
        𝜃ₙ;𝑪;𝜺;𝜿;𝑵;𝑴;
    )
end

function generateMsh(filepath::String; lc = 1.0, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Patch Test")

    Ω = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    tag = BenchmarkExample.addEdgeElements((2,1), order)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(𝑅*𝜃, 0.0, 0.0, lc, 2)
    𝐴 = gmsh.model.geo.addPoint(𝑅*𝜃,  𝐿, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0, 𝐿, 0.0, lc, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [1,2,3,4], -1, "Γ")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Ω
end
end