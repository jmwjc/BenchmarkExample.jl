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

end
