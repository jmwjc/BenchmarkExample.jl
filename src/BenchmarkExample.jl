module BenchmarkExample

import Gmsh: gmsh
include("ScordelisLoRoof.jl")

function addEdgeElements(dimTag::Tuple{Int,Int}, order::Int=1)
    dim, tag = dimTag
    gmsh.model.mesh.createEdges(dimTag)
    elementTypes, ~, nodeTags = gmsh.model.mesh.getElements(dim,tag)
    s = gmsh.model.addDiscreteEntity(dim-1)
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType,tag)
        edgeTags, edgeOrientations = gmsh.model.mesh.getEdges(edgeNodes)
        maxTag = gmsh.model.mesh.getMaxElementTag()
        elementTags = [i+maxTag for i in 1:length(edgeTags)]
        type = gmsh.model.mesh.getElementType("Line", order)
        gmsh.model.mesh.addElementsByType(s, type, elementTags, edgeNodes)
    end
    return (dim-1,s)
end

function cylindricalCoordinate(𝑅::Float64)
    𝒂₁(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = cos(ξ¹/𝑅), 0.0, -sin(ξ¹/𝑅)
    𝒂₂(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0, 1.0, 0.0
    𝒂₃(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = sin(ξ¹/𝑅), 1.0, cos(ξ¹/𝑅)
    a¹¹(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 1.0
    a²²(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 1.0
    a³³(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 1.0
    a¹²(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    a¹³(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    a²³(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    Γ¹₁₁(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    Γ²₁₁(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    Γ¹₁₂(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    Γ²₁₂(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    Γ¹₂₂(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    Γ²₂₂(ξ¹::Float64,ξ²::Float64,ξ³::Float64) = 0.0
    return ()->(𝒂₁;𝒂₂;𝒂₃;a¹¹;a²²;a³³;a¹²;a¹³;a²³;Γ¹₁₁;Γ²₁₁;Γ¹₁₂;Γ²₁₂;Γ¹₂₂;Γ²₂₂)
end

end
