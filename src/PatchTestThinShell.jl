module PatchTestThinShell

import ..BenchmarkExample
import Gmsh: gmsh
import Tensors: ⋅, ⊗, ×, ⊡, Vec, Tensor, SymmetricTensor, gradient, divergence, curl

𝐸 = 1.0
𝜈 = 0.0
𝑅 = 1.0
# 𝜃 = 0.1*π
𝜃 = 1.0
𝐿 = 𝑅*𝜃
ℎ = 0.05
# ℎ = 12.0^0.5

function variables(cs::Function,𝒖::Function)
    function ∇𝒖(x)
        x_ = Vec{3}((x[1],x[2],0.0))
        gradient(𝒖,x_)
    end
    function ∇²𝒖(x)
        x_ = Vec{3}((x[1],x[2],0.0))
        gradient(∇𝒖,x_)
    end
    function 𝜺(x)
        grad𝒖 = ∇𝒖(x)
        ∂₁𝒖 = Tensor{1,3}((grad𝒖[:,1]))
        ∂₂𝒖 = Tensor{1,3}((grad𝒖[:,2]))
        ε₁₁ = cs.𝒂₁(x)⋅∂₁𝒖
        ε₂₂ = cs.𝒂₂(x)⋅∂₂𝒖
        ε₁₂ = 0.5*(cs.𝒂₁(x)⋅∂₂𝒖 + cs.𝒂₂(x)⋅∂₁𝒖)
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
        Γ¹₂₂ = cs.Γ¹₂₂(x)
        Γ²₁₁ = cs.Γ²₁₁(x)
        Γ²₁₂ = cs.Γ²₁₂(x)
        Γ²₂₂ = cs.Γ²₂₂(x)
        𝒂₃ = cs.𝒂₃(x)
        κ₁₁ = (Γ¹₁₁.*∂₁𝒖+Γ²₁₁.*∂₂𝒖 - ∂₁₁𝒖)⋅𝒂₃
        κ₂₂ = (Γ¹₂₂.*∂₁𝒖+Γ²₂₂.*∂₂𝒖 - ∂₂₂𝒖)⋅𝒂₃
        κ₁₂ = (Γ¹₁₂.*∂₁𝒖+Γ²₁₂.*∂₂𝒖 - ∂₁₂𝒖)⋅𝒂₃
        return Tensor{2,2}((κ₁₁,κ₁₂,κ₁₂,κ₂₂))
    end
    function 𝑪(x)
        a¹¹ = cs.a¹¹(x)
        a²² = cs.a²²(x)
        a¹² = cs.a¹²(x)
        C¹¹¹¹ = a¹¹*a¹¹
        C²²²² = a²²*a²²
        C¹¹²² = 𝜈*a¹¹*a²² + (1-𝜈)*a¹²*a¹²
        C¹¹¹² = a¹¹*a¹²
        C²²¹² = a²²*a¹²
        C¹²¹² = 0.5*((1-𝜈)*a¹¹*a²² + (1+𝜈)*a¹²*a¹²)
        return SymmetricTensor{4,2}((C¹¹¹¹,C¹¹¹²,C¹¹²²,C¹¹¹²,C¹²¹²,C²²¹²,C¹¹²²,C²²¹²,C²²²²))
    end
    function 𝑵(x)
        Dᵐ = 𝐸*ℎ
        ε = 𝜺(x)
        C = 𝑪(x)
        return Tensor{2,2}((Dᵐ.*(C ⊡ ε)))
    end
    function 𝑴(x)
        Dᵇ = 𝐸*ℎ^3/12
        κ = 𝜿(x)
        C = 𝑪(x)
        return Tensor{2,2}((Dᵇ.*(C ⊡ κ)))
    end
    function 𝜃ₙ(x,n¹,n²)
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
        function ∇𝒂₃_(x_)
            x__ = Vec{3}((x_[1],x_[2],0.0))
            return gradient(cs.𝒂₃,x__)
        end
        ∇𝒂₃ = ∇𝒂₃_(x)
        x__ = Vec{3}((x[1],x[2],0.0))
        ∇²𝒂₃ = gradient(∇𝒂₃_,x__)
        ∂₁𝒂₃ = ∇𝒂₃_(x)[:,1]
        ∂₂𝒂₃ = ∇𝒂₃_(x)[:,2]
        Γ¹₁₁ = cs.Γ¹₁₁(x)
        Γ¹₁₂ = cs.Γ¹₁₂(x)
        Γ¹₂₂ = cs.Γ¹₂₂(x)
        Γ²₁₁ = cs.Γ²₁₁(x)
        Γ²₂₁ = cs.Γ²₁₂(x)
        Γ²₂₂ = cs.Γ²₂₂(x)
        ∂₁Γ¹₁₁ = cs.∂₁Γ¹₁₁(x)
        ∂₂Γ¹₁₁ = cs.∂₂Γ¹₁₁(x)
        ∂₁Γ¹₁₂ = cs.∂₁Γ¹₁₂(x)
        ∂₂Γ¹₁₂ = cs.∂₂Γ¹₁₂(x)
        ∂₁Γ¹₂₂ = cs.∂₁Γ¹₂₂(x)
        ∂₂Γ¹₂₂ = cs.∂₂Γ¹₂₂(x)
        ∂₁Γ²₁₁ = cs.∂₁Γ²₁₁(x)
        ∂₂Γ²₁₁ = cs.∂₂Γ²₁₁(x)
        ∂₁Γ²₁₂ = cs.∂₁Γ²₁₂(x)
        ∂₂Γ²₁₂ = cs.∂₂Γ²₁₂(x)
        ∂₁Γ²₂₂ = cs.∂₁Γ²₂₂(x)
        ∂₂Γ²₂₂ = cs.∂₂Γ²₂₂(x)
        Γ¹ = Tensor{2,2}((Γ¹₁₁,Γ¹₁₂,Γ¹₁₂,Γ¹₂₂))
        Γ² = Tensor{2,2}((Γ²₁₁,Γ²₂₁,Γ²₂₁,Γ²₂₂))
        b = Tensor{2,2}((b₁₁,b₁₂,b₁₂,b₂₂))
        Γᵞᵧ = Vec{2}((Γ¹₁₁+Γ²₂₁,Γ¹₁₂+Γ²₂₂))
        ∂ᵧΓᵞ = Tensor{2,2}((∂₁Γ¹₁₁+∂₂Γ²₁₁,∂₁Γ¹₁₂+∂₂Γ²₁₂,∂₁Γ¹₁₂+∂₂Γ²₁₂,∂₁Γ¹₂₂+∂₂Γ²₂₂))
        ∂Γᵞᵧ = Tensor{2,2}((∂₁Γ¹₁₁+∂₁Γ²₁₂,∂₂Γ¹₁₁+∂₂Γ²₁₂,∂₁Γ¹₁₂+∂₁Γ²₂₂,∂₂Γ¹₁₂+∂₂Γ²₂₂))
        N = 𝑵(x)
        ∇N = gradient(𝑵,x)
        divN = ∇N[:,1,1] + ∇N[:,2,2]
        NΓᵞᵧ = N ⋅ Γᵞᵧ
        M = 𝑴(x)
        ∇𝑴(x_) = gradient(𝑴,x_)
        ∇M = ∇𝑴(x)
        ∂₁𝑴(x_) = Tensor{2,2}(∇𝑴(x_)[:,:,1])
        ∂₂𝑴(x_) = Tensor{2,2}(∇𝑴(x_)[:,:,2])
        ∂₁M = ∂₁𝑴(x)
        ∂₂M = ∂₂𝑴(x)
        ∇∂₁𝑴(x_) = gradient(∂₁𝑴,x_)
        ∇∂₂𝑴(x_) = gradient(∂₂𝑴,x_)
        ∇∂₁M = ∇∂₁𝑴(x)
        ∇∂₂M = ∇∂₂𝑴(x)
        divM = ∇M[:,1,1] + ∇M[:,2,2]
        div2M = ∇∂₁M[1,1,1]+∇∂₂M[2,1,1]+∇∂₁M[1,2,2]+∇∂₂M[2,2,2]
        ∇²𝒂₃M = ∇²𝒂₃[:,1,1]*M[1,1]+∇²𝒂₃[:,1,2]*M[1,2]+∇²𝒂₃[:,2,1]*M[2,1]+∇²𝒂₃[:,2,2]*M[2,2]
        ΓᵞᵧM = Γᵞᵧ⋅M
        𝒃1 = (Γ¹ ⊡ N + divN[1] + NΓᵞᵧ[1]).*𝒂₁ + (Γ² ⊡ N + divN[2] + NΓᵞᵧ[2]).*𝒂₂ + b ⊡ N.*𝒂₃
        𝒃2 = (∂ᵧΓᵞ ⊡ M .*𝒂₃ + Γ¹ ⊡ M .* ∂₁𝒂₃ + Γ¹ ⊡ ∂₁M .* 𝒂₃ + Γ² ⊡ M .* ∂₂𝒂₃ + Γ² ⊡ ∂₂M .* 𝒂₃
           + Γᵞᵧ[1]*(Γ¹ ⊡ M) .* 𝒂₃ + Γᵞᵧ[2]*(Γ² ⊡ M) .* 𝒂₃ + div2M.*𝒂₃ + 2*(∇𝒂₃[:,1].*divM[1] + ∇𝒂₃[:,2].*divM[2]) + ∇²𝒂₃M + ∂Γᵞᵧ ⊡ M .* 𝒂₃
           + 2*(∇𝒂₃[:,1].*ΓᵞᵧM[1]+∇𝒂₃[:,2].*ΓᵞᵧM[2] + Γᵞᵧ⋅divM.*𝒂₃) + Γᵞᵧ⋅M⋅Γᵞᵧ.*𝒂₃)
        return -𝒃1-𝒃2
        # return -𝒃1
        # return -𝒃2
    end
    return ()->(
        𝜃ₙ;𝑪;𝜺;𝜿;𝑵;𝑴;𝒃;∇𝒖;
    )
end

function generateMsh(filepath::String; lc = 1.0, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Patch Test")

    Ω = generateGeo(lc)

    if quad
        gmsh.model.mesh.setRecombine(2, Ω)
    end
    
    gmsh.model.mesh.setAlgorithm(2, Ω, 1)
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
    gmsh.model.geo.addPoint(𝑅*𝜃,  𝐿, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0, 𝐿, 0.0, lc, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    Ω = gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [1], -1, "Γ¹")
    gmsh.model.addPhysicalGroup(1, [2], -1, "Γ²")
    gmsh.model.addPhysicalGroup(1, [3], -1, "Γ³")
    gmsh.model.addPhysicalGroup(1, [4], -1, "Γ⁴")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    return Ω
end
end