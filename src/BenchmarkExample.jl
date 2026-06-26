module BenchmarkExample

import Gmsh: gmsh
import Tensors: ⋅, ⊗, ×, Vec, gradient, divergence, curl
include("T3.jl")
include("PatchTest.jl")
include("CookMembrance.jl")
include("CantileverBeam.jl")
include("PatchTestThinShell.jl")
include("ScordelisLoRoof.jl")
include("SphericalShell.jl")
include("SquarePlate.jl")
include("Circular.jl")
include("MorleysAcuteSkewPlate.jl")
include("PlateWithHole.jl")
include("Punch.jl")
include("Block_plot.jl")
include("Block.jl")
include("PatchTest3D.jl")
include("Tet4.jl")
include("TickWalledCylinder.jl")
include("SquarePrism.jl")
include("Balloon.jl")
include("Balloon2.jl")
include("Waterstop.jl")
include("Isolation.jl")
include("Isolation_half.jl")
# function addEdgeElements(dimTag::Tuple{Int,Int}, order::Int=1)
#     dim, tag = dimTag
#     gmsh.model.mesh.createEdges(dimTag)
#     elementTypes, ~, nodeTags = gmsh.model.mesh.getElements(dim,tag)
#     s = gmsh.model.addDiscreteEntity(1)
    
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)
#         edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true) #nodetags
#         edgeTags, edgeOrientations = gmsh.model.mesh.getEdges(edgeNodes)      #elementtags
#         maxTag = gmsh.model.mesh.getMaxElementTag()
#         elementTags = [i+maxTag+length(edgeNodes) for i in 1:length(edgeTags)]
#         type = gmsh.model.mesh.getElementType("Line", order)
#         edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType,tag)
#         gmsh.model.mesh.addElementsByType(s, type, elementTags, edgeNodes)
#     end
  
  
#     return s
# end

function addEdgeElements(dimTag::Tuple{Int,Int}, order::Int=1)
    dim, tag = dimTag
    # gmsh.model.mesh.createEdges(dimTag)
    elementTypes, ~, nodeTags = gmsh.model.mesh.getElements(dim,tag)
    s = gmsh.model.addDiscreteEntity(1)
    
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true) #nodetags
        # edgeTags, edgeOrientations = gmsh.model.mesh.getEdges(edgeNodes)      #elementtags
        
        # elementTags = [i+maxTag+length(edgeNodes) for i in 1:length(edgeTags)]
        type = gmsh.model.mesh.getElementType("Line", order)
        
        gmsh.model.mesh.addElementsByType(s, type, Int[], edgeNodes)
    end
  
  
    return s
end


# function addFaceElements(dimTag::Tuple{Int,Int}, order::Int=1)
#     dim, tag = dimTag
#     gmsh.model.mesh.createFaces(dimTag)
#     elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim,tag)
#     s = gmsh.model.addDiscreteEntity(2)
#     faceTypeMap = Dict(
#         4 => (faceType=3, faceName="Triangle"),  # 四面体 → 三角形
#         5 => (faceType=4, faceName="Quadrangle") # 六面体 → 四边形
#     )
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)

#         faceType = faceTypeMap[elementType].faceType
#         faceName = faceTypeMap[elementType].faceName
#         faceNodes = gmsh.model.mesh.getElementFaceNodes(elementType,faceType,tag)
#         faceTags, faceOrientations = gmsh.model.mesh.getFaces(faceType,faceNodes)
#         maxTag = gmsh.model.mesh.getMaxElementTag()
#         newTags = collect(maxTag+1 : maxTag+length(faceTags))
#         # newTags = [i+maxTag+length(faceNodes) for i in 1:length(faceTags)]
#         type = gmsh.model.mesh.getElementType(faceName, order)
#         gmsh.model.mesh.addElementsByType(s, type, newTags, faceNodes)
#     end
#     return s
# end

# function addFaceElements(dimTag::Tuple{Int,Int}, order::Int=1)
#     dim, tag = dimTag
#     # gmsh.model.mesh.createFaces(dimTag)
#     # 1. 获取体积单元信息
#     elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim,tag)
#     s = gmsh.model.addDiscreteEntity(2)
#     faceTypeMap = Dict(
#         4 => (faceType=3, faceName="Triangle",nf = 3),  # 四面体 → 三角形
#         5 => (faceType=4, faceName="Quadrangle",nf = 4) # 六面体 → 四边形
#     )
#     maxTag = gmsh.model.mesh.getMaxElementTag()
    
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)

#         faceType = faceTypeMap[elementType].faceType
#         faceName = faceTypeMap[elementType].faceName
#         faceNodes = gmsh.model.mesh.getElementFaceNodes(elementType,faceType,tag)
#         # faceTags, faceOrientations = gmsh.model.mesh.getFaces(faceType,faceNodes)
        
#         type = gmsh.model.mesh.getElementType(faceName, order)
#         gmsh.model.mesh.addElementsByType(s, type, Int[], faceNodes)
#     end
#     return s
# end


function addFaceElements(dimTag::Tuple{Int,Int}, order::Int=1)
    dim, tag = dimTag
    @assert dim == 3 "addFaceElements only supports 3D volume entities"

    # volume elements on this entity
    elementTypes, elementTagsVec, nodeTagsVec = gmsh.model.mesh.getElements(dim, tag)

    # create one discrete surface entity to hold all local faces
    s = gmsh.model.addDiscreteEntity(2)

    # only for duplicated local faces, no deduplication
    faceInfo = Dict(
        4 => (faceType = 3, faceName = "Triangle",    nb = 4, nf = 3), # Tet4
        5 => (faceType = 4, faceName = "Quadrangle",  nb = 6, nf = 4)  # Hex8
    )

    for elementType in elementTypes
        @assert haskey(faceInfo, elementType) "Unsupported volume element type: $elementType"

        info = faceInfo[elementType]

        # Gmsh doc meaning:
        # elementType = volume type code
        # faceType    = number of nodes on each face (3 for tri face, 4 for quad face)
        faceNodes = gmsh.model.mesh.getElementFaceNodes(elementType, info.faceType, tag)

        faceElementType = gmsh.model.mesh.getElementType(info.faceName, order)

        # Keep duplicated faces intentionally:
        # every volume element contributes all its local faces
        gmsh.model.mesh.addElementsByType(s, faceElementType, Int[], faceNodes)
    end

    return s
end

function debug_face_parent_mapping(dimTag::Tuple{Int,Int}, surfaceTag::Int)
    dim, volTag = dimTag
    @assert dim == 3 "debug_face_parent_mapping only supports 3D volume entities"

    # -----------------------------
    # 1) volume elements on this volume entity
    # -----------------------------
    elementTypesΩ, elementTagsΩ_vec, nodeTagsΩ_vec = gmsh.model.mesh.getElements(dim, volTag)

    println("===================================================")
    println("Volume entity tag = ", volTag)
    println("Surface entity tag = ", surfaceTag)
    println("===================================================")

    # Build a dictionary:
    #   key = sorted face node tuple
    #   value = (parent volume element tag, local face id, volume element type)
    face_to_parent = Dict{Tuple{Vararg{UInt64}}, Tuple{UInt64,Int,Int32}}()

    for (elementTypeΩ, elementTagsΩ, nodeTagsΩ) in zip(elementTypesΩ, elementTagsΩ_vec, nodeTagsΩ_vec)
        _, _, _, niΩ = gmsh.model.mesh.getElementProperties(elementTypeΩ)
        neΩ = Int(length(nodeTagsΩ) ÷ niΩ)

        # number of nodes per face argument for getElementFaceNodes
        faceTypeArg =
            elementTypeΩ == 4 ? 3 :   # Tet4 faces are Tri3
            elementTypeΩ == 5 ? 4 :   # Hex8 faces are Quad4
            error("Unsupported volume elementTypeΩ = $elementTypeΩ")

        nb =
            elementTypeΩ == 4 ? 4 :
            elementTypeΩ == 5 ? 6 :
            error("Unsupported volume elementTypeΩ = $elementTypeΩ")

        allFaceNodes = gmsh.model.mesh.getElementFaceNodes(elementTypeΩ, faceTypeArg, volTag)

        # allFaceNodes is grouped by volume element, then by local face
        for e in 1:neΩ
            tagΩ = elementTagsΩ[e]

            for lf in 1:nb
                start = faceTypeArg * (nb*(e-1) + (lf-1)) + 1
                stop  = start + faceTypeArg - 1
                face_nodes = allFaceNodes[start:stop]
                key = Tuple(sort(collect(face_nodes)))

                face_to_parent[key] = (tagΩ, lf, elementTypeΩ)
            end
        end
    end

    # -----------------------------
    # 2) read generated face elements on the discrete surface entity
    # -----------------------------
    elementTypesS, elementTagsS_vec, nodeTagsS_vec = gmsh.model.mesh.getElements(2, surfaceTag)

    println("Generated face elements on surface entity:")
    total_faces = 0

    for (elementTypeS, elementTagsS, nodeTagsS) in zip(elementTypesS, elementTagsS_vec, nodeTagsS_vec)
        _, _, _, niS = gmsh.model.mesh.getElementProperties(elementTypeS)
        neS = Int(length(nodeTagsS) ÷ niS)

        println("  surface elementType = ", elementTypeS, ", niS = ", niS, ", neS = ", neS)

        for k in 1:neS
            total_faces += 1
            tagS = elementTagsS[k]
            face_nodes = nodeTagsS[niS*(k-1)+1 : niS*k]
            key = Tuple(sort(collect(face_nodes)))

            if haskey(face_to_parent, key)
                parentTag, localFaceId, parentType = face_to_parent[key]
                println(
                    "k = ", total_faces,
                    " | surfaceTag = ", tagS,
                    " | nodes = ", face_nodes,
                    " | parent volume tag = ", parentTag,
                    " | local face id = ", localFaceId,
                    " | parent elementType = ", parentType
                )
            else
                println(
                    "k = ", total_faces,
                    " | surfaceTag = ", tagS,
                    " | nodes = ", face_nodes,
                    " | parent = NOT FOUND"
                )
            end
        end
    end

    println("===================================================")
    println("Total generated faces checked = ", total_faces)
    println("===================================================")
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
    function 𝚪₁₁(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(𝒂₁,𝛏_)[:,1]
    end
    function 𝚪₂₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(𝒂₂,𝛏_)[:,2]
    end
    function 𝚪₁₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(𝒂₁,𝛏_)[:,2]
    end
    b₁₁(𝛏::Vec) = 𝒂₃(𝛏)⋅𝚪₁₁(𝛏)
    b₂₂(𝛏::Vec) = 𝒂₃(𝛏)⋅𝚪₂₂(𝛏)
    b₁₂(𝛏::Vec) = 𝒂₃(𝛏)⋅𝚪₁₂(𝛏)
    Γ¹₁₁(𝛏::Vec) = 𝒂¹(𝛏)⋅𝚪₁₁(𝛏)
    Γ²₁₁(𝛏::Vec) = 𝒂²(𝛏)⋅𝚪₁₁(𝛏)
    Γ¹₂₂(𝛏::Vec) = 𝒂¹(𝛏)⋅𝚪₂₂(𝛏)
    Γ²₂₂(𝛏::Vec) = 𝒂²(𝛏)⋅𝚪₂₂(𝛏)
    Γ¹₁₂(𝛏::Vec) = 𝒂¹(𝛏)⋅𝚪₁₂(𝛏)
    Γ²₁₂(𝛏::Vec) = 𝒂²(𝛏)⋅𝚪₁₂(𝛏)
    function ∂₁Γ¹₁₁(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ¹₁₁,𝛏_)[1]
    end
    function ∂₂Γ¹₁₁(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ¹₁₁,𝛏_)[2]
    end
    function ∂₁Γ¹₂₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ¹₂₂,𝛏_)[1]
    end
    function ∂₂Γ¹₂₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ¹₂₂,𝛏_)[2]
    end
    function ∂₁Γ¹₁₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ¹₁₂,𝛏_)[1]
    end
    function ∂₂Γ¹₁₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ¹₁₂,𝛏_)[2]
    end
    function ∂₁Γ²₁₁(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ²₁₁,𝛏_)[1]
    end
    function ∂₂Γ²₁₁(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ²₁₁,𝛏_)[2]
    end
    function ∂₁Γ²₂₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ²₂₂,𝛏_)[1]
    end
    function ∂₂Γ²₂₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ²₂₂,𝛏_)[2]
    end
    function ∂₁Γ²₁₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ²₁₂,𝛏_)[1]
    end
    function ∂₂Γ²₁₂(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(Γ²₁₂,𝛏_)[2]
    end
    function ∂₁𝒂₃(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(𝒂₃,𝛏_)[:,1]
    end
    function ∂₂𝒂₃(𝛏::Vec)
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(𝒂₃,𝛏_)[:,2]
    end
    function ∂₁₁𝒂₃(𝛏::Vec) 
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(∂₁𝒂₃,𝛏_)[:,1]
    end
    function ∂₁₂𝒂₃(𝛏::Vec) 
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(∂₁𝒂₃,𝛏_)[:,2]
    end
    function ∂₂₂𝒂₃(𝛏::Vec) 
        𝛏_ = Vec{3}((𝛏[1],𝛏[2],0.0))
        return gradient(∂₂𝒂₃,𝛏_)[:,2]
    end
    𝐽(𝛏::Vec) = 𝒂₁(𝛏)⋅(𝒂₂(𝛏)×𝒂₃(𝛏))
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
        ∂₁₁𝒂₃;∂₂₂𝒂₃;∂₁₂𝒂₃;𝐽;
    )
end
function sphericalCoordinate(𝑅::Float64)
    𝒂₁(𝛏::Vec) = Vec{3}((-sin(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅),  cos(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅), 0.0))
    𝒂₂(𝛏::Vec) = Vec{3}((-cos(𝛏[1]/𝑅)*sin(𝛏[2]/𝑅), -sin(𝛏[1]/𝑅)*sin(𝛏[2]/𝑅), cos(𝛏[2]/𝑅)))
    𝒂₃(𝛏::Vec) = Vec{3}(( cos(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅),  sin(𝛏[1]/𝑅)*cos(𝛏[2]/𝑅), sin(𝛏[2]/𝑅)))
    𝒂¹(𝛏::Vec) = Vec{3}((-sin(𝛏[1]/𝑅)*sec(𝛏[2]/𝑅),  cos(𝛏[1]/𝑅)*sec(𝛏[2]/𝑅), 0.0))
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
    ∂₁₁𝒂₃;∂₂₂𝒂₃;∂₁₂𝒂₃;𝐽;
)
end
function cartesianCoordinate()
    𝒂₁(𝛏::Vec) = Vec{3}((𝛏[1]^0, 0.0, 0.0))
    𝒂₂(𝛏::Vec) = Vec{3}((0.0, 𝛏[1]^0, 0.0))
    𝒂₃(𝛏::Vec) = Vec{3}((0.0, 0.0, 𝛏[1]^0))
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
        ∂₁₁𝒂₃;∂₂₂𝒂₃;∂₁₂𝒂₃;𝐽;
    )
end
end
end
