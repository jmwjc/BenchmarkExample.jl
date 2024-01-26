module BenchmarkExample

import Gmsh: gmsh
import Tensors: ⋅, ⊗, Vec, gradient, divergence, curl

include("ScordelisLoRoof.jl")
include("SphericalShell.jl")

function addEdgeElements(dimTag::Tuple{Int,Int}, order::Int=1)
    dim, tag = dimTag
    gmsh.model.mesh.createEdges(dimTag)
    elementTypes, ~, nodeTags = gmsh.model.mesh.getElements(dim,tag)
    s = gmsh.model.addDiscreteEntity(1)
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType,tag)
        edgeTags, edgeOrientations = gmsh.model.mesh.getEdges(edgeNodes)
        maxTag = gmsh.model.mesh.getMaxElementTag()
        elementTags = [i+maxTag+length(edgeNodes) for i in 1:length(edgeTags)]
        type = gmsh.model.mesh.getElementType("Line", order)
        gmsh.model.mesh.addElementsByType(s, type, elementTags, edgeNodes)
    end
    return s
end



const ∇ = Val(:∇)
⊗(∇::Val{:∇},f::Function) = (x)->gradient(f,x)
⋅(∇::Val{:∇},f::Function) = (x)->divergence(f,x)
×(∇::Val{:∇},f::Function) = (x)->curl(f,x)

covariantDerivative = quote
    a¹¹(𝛏::Vec) = 𝒂¹(𝛏)⋅𝒂¹(𝛏)
    a¹²(𝛏::Vec) = 𝒂¹(𝛏)⋅𝒂²(𝛏)
    a¹³(𝛏::Vec) = 𝒂¹(𝛏)⋅𝒂³(𝛏)
    a²²(𝛏::Vec) = 𝒂²(𝛏)⋅𝒂²(𝛏)
    a²³(𝛏::Vec) = 𝒂²(𝛏)⋅𝒂³(𝛏)
    a³³(𝛏::Vec) = 𝒂³(𝛏)⋅𝒂³(𝛏)
    a₁₁(𝛏::Vec) = 𝒂₁(𝛏)⋅𝒂₁(𝛏)
    a₁₂(𝛏::Vec) = 𝒂₁(𝛏)⋅𝒂₂(𝛏)
    a₁₃(𝛏::Vec) = 𝒂₁(𝛏)⋅𝒂₃(𝛏)
    a₂₂(𝛏::Vec) = 𝒂₂(𝛏)⋅𝒂₂(𝛏)
    a₂₃(𝛏::Vec) = 𝒂₂(𝛏)⋅𝒂₃(𝛏)
    a₃₃(𝛏::Vec) = 𝒂₃(𝛏)⋅𝒂₃(𝛏)
    𝚪₁₁(𝛏::Vec) = gradient(𝒂₁,𝛏)[:,1]
    𝚪₂₂(𝛏::Vec) = gradient(𝒂₂,𝛏)[:,2]
    𝚪₁₂(𝛏::Vec) = gradient(𝒂₁,𝛏)[:,2]
    b₁₁(𝛏::Vec) = 𝒂₃(𝛏)⋅𝚪₁₁(𝛏)
    b₂₂(𝛏::Vec) = 𝒂₃(𝛏)⋅𝚪₂₂(𝛏)
    b₁₂(𝛏::Vec) = 𝒂₃(𝛏)⋅𝚪₁₂(𝛏)
    Γ¹₁₁(𝛏::Vec) = 𝒂¹(𝛏)⋅𝚪₁₁(𝛏)
    Γ²₁₁(𝛏::Vec) = 𝒂²(𝛏)⋅𝚪₁₁(𝛏)
    Γ¹₂₂(𝛏::Vec) = 𝒂¹(𝛏)⋅𝚪₂₂(𝛏)
    Γ²₂₂(𝛏::Vec) = 𝒂²(𝛏)⋅𝚪₂₂(𝛏)
    Γ¹₁₂(𝛏::Vec) = 𝒂¹(𝛏)⋅𝚪₁₂(𝛏)
    Γ²₁₂(𝛏::Vec) = 𝒂²(𝛏)⋅𝚪₁₂(𝛏)
    ∂₁Γ¹₁₁(𝛏::Vec) = gradient(Γ¹₁₁,𝛏)[1]
    ∂₂Γ¹₁₁(𝛏::Vec) = gradient(Γ¹₁₁,𝛏)[2]
    ∂₁Γ¹₂₂(𝛏::Vec) = gradient(Γ¹₂₂,𝛏)[1]
    ∂₂Γ¹₂₂(𝛏::Vec) = gradient(Γ¹₂₂,𝛏)[2]
    ∂₁Γ¹₁₂(𝛏::Vec) = gradient(Γ¹₁₂,𝛏)[1]
    ∂₂Γ¹₁₂(𝛏::Vec) = gradient(Γ¹₁₂,𝛏)[2]
    ∂₁Γ²₁₁(𝛏::Vec) = gradient(Γ²₁₁,𝛏)[1]
    ∂₂Γ²₁₁(𝛏::Vec) = gradient(Γ²₁₁,𝛏)[2]
    ∂₁Γ²₂₂(𝛏::Vec) = gradient(Γ²₂₂,𝛏)[1]
    ∂₂Γ²₂₂(𝛏::Vec) = gradient(Γ²₂₂,𝛏)[2]
    ∂₁Γ²₁₂(𝛏::Vec) = gradient(Γ²₁₂,𝛏)[1]
    ∂₂Γ²₁₂(𝛏::Vec) = gradient(Γ²₁₂,𝛏)[2]
    ∂₁𝒂₃(𝛏::Vec) = gradient(𝒂₃,𝛏)[:,1]
    ∂₂𝒂₃(𝛏::Vec) = gradient(𝒂₃,𝛏)[:,2]
    ∂₁₁𝒂₃(𝛏::Vec) = gradient(∂₁𝒂₃,𝛏)[:,1]
    ∂₁₂𝒂₃(𝛏::Vec) = gradient(∂₁𝒂₃,𝛏)[:,2]
    ∂₂₂𝒂₃(𝛏::Vec) = gradient(∂₂𝒂₃,𝛏)[:,2]
end

@eval begin
function cylindricalCoordinate(𝑅::Float64)
    𝒂₁(𝛏::Vec) = Vec{3}((cos(𝛏[1]/𝑅), 0.0, -sin(𝛏[1]/𝑅)))
    𝒂₂(𝛏::Vec) = Vec{3}((0.0, 1.0, 0.0))
    𝒂₃(𝛏::Vec) = Vec{3}((sin(𝛏[1]/𝑅), 0.0, cos(𝛏[1]/𝑅)))
    𝒂¹(𝛏::Vec) = 𝒂₁(𝛏)
    𝒂²(𝛏::Vec) = 𝒂₂(𝛏)
    𝒂³(𝛏::Vec) = 𝒂₃(𝛏)

    $covariantDerivative

    return ()->(
        𝒂₁;𝒂₂;𝒂₃;
        𝒂¹;𝒂²;𝒂³;
        a¹¹;a²²;a¹²;
        a₁₁;a₂₂;a₁₂;
        b₁₁;b₂₂;b₁₂;
        Γ¹₁₁;Γ²₁₁;Γ¹₂₂;Γ²₂₂;Γ¹₁₂;Γ²₁₂;
        ∂₁Γ¹₁₁;∂₁Γ²₁₁;∂₁Γ¹₂₂;∂₁Γ²₂₂;∂₁Γ¹₁₂;∂₁Γ²₁₂;
        ∂₂Γ¹₁₁;∂₂Γ²₁₁;∂₂Γ¹₂₂;∂₂Γ²₂₂;∂₂Γ¹₁₂;∂₂Γ²₁₂;
        ∂₁₁𝒂₃;∂₂₂𝒂₃;∂₁₂𝒂₃
    )
end
function sphericalCoordinate(𝑅::Float64)
    𝒂₁(𝛏::Vec) = Vec{3}((-sin(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅), cos(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅), 0.0))
    𝒂₂(𝛏::Vec) = Vec{3}((-cos(𝛏[1]/𝑅)*sin(𝛏[2]/𝑅), -sin(𝛏[1]/𝑅)*sin(𝛏[2]/𝑅), cos(𝛏[2]/𝑅)))
    𝒂₃(𝛏::Vec) = Vec{3}((cos(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅), sin(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅), sin(𝛏[2]/𝑅)))
    𝒂¹(𝛏::Vec) = Vec{3}((-sin(𝛏[1]/𝑅)*sec(𝛏[2]/𝑅), cos(𝛏[1]/𝑅)*sec(𝛏[2]/𝑅), 0.0))
    𝒂²(𝛏::Vec) =  𝒂₂(𝛏)
    𝒂³(𝛏::Vec) =  𝒂₃(𝛏)

    $covariantDerivative

    return ()->(
    𝒂₁;𝒂₂;𝒂₃;
    𝒂¹;𝒂²;𝒂³;
    a¹¹;a²²;a¹²;
    a₁₁;a₂₂;a₁₂;
    b₁₁;b₂₂;b₁₂;
    Γ¹₁₁;Γ²₁₁;Γ¹₂₂;Γ²₂₂;Γ¹₁₂;Γ²₁₂;
    ∂₁Γ¹₁₁;∂₁Γ²₁₁;∂₁Γ¹₂₂;∂₁Γ²₂₂;∂₁Γ¹₁₂;∂₁Γ²₁₂;
    ∂₂Γ¹₁₁;∂₂Γ²₁₁;∂₂Γ¹₂₂;∂₂Γ²₂₂;∂₂Γ¹₁₂;∂₂Γ²₁₂;
    ∂₁₁𝒂₃;∂₂₂𝒂₃;∂₁₂𝒂₃
)
end
end
end
