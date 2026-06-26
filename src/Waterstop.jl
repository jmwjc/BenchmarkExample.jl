module Waterstop

import ..BenchmarkExample
import Gmsh: gmsh

# --- 1. 几何参数 (5腿4拱, 侧面外凸) ---
const H_total = 22.0
const W_top_half = 20.5

# 底部结构
const N_legs = 5            # 5条腿
const W_leg = 6.0           # 腿宽
const W_arch = 6.0          # 拱门宽 (间隙)
const H_arch_straight = 4.5 # 拱门直立高度

# 侧面飞翼
const H_wing_start = 7.5    # 腿部高度
const H_wing_tip = 9.5      # 飞翼尖端高度
const W_wing_stickout = 3.5 # 侧面外凸程度

# 孔洞
const R_hole_large = 3.2
const R_hole_small = 2.5
const Y_hole_top = 16.5
const Y_hole_bot = 11.5

function generateMsh(filepath::String; lc = 0.8, transfinite = -1, order = 1, quad = false)
    gmsh.initialize()
    gmsh.model.add("Waterstop")

    # 1. 生成几何
    S1, boundary_lines = generateGeo(lc)

    # 【修复重点 1】: 必须先同步几何，才能进行网格划分或添加物理组
    # gmsh.model.geo.synchronize()

    # 2. 网格设置
    if quad
        gmsh.model.mesh.setRecombine(2, S1) 
    end
    
    
    gmsh.model.mesh.setAlgorithm(2, S1, 5)
     gmsh.model.mesh.generate(2)
    # 【修复重点 2】: 使用 S1 变量，而不是硬编码的 1
    # 之前的 (2,1) 会导致越界，因为 S1 可能不等于 1
    gmsh.model.mesh.SecondOrderLinear= 1
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
    gmsh.model.mesh.setOrder(order)
 gmsh.model.mesh.refine()


    tag = BenchmarkExample.addEdgeElements((2, 1), order)

    # 为 addEdgeElements 返回的边界添加物理组
    gmsh.model.addPhysicalGroup(1, [tag], -1, "Γ")

    gmsh.model.geo.synchronize()
    
    # 最后写入
    gmsh.write(filepath)
    gmsh.finalize()
end

@inline function generateGeo(lc = 1.0)
    
    # --- 辅助函数：画圆孔 ---
    function addHole(x, y, r, lc)
        p_cen = gmsh.model.geo.addPoint(x, y, 0, lc)
        p1 = gmsh.model.geo.addPoint(x + r, y, 0, lc)
        p2 = gmsh.model.geo.addPoint(x, y + r, 0, lc)
        p3 = gmsh.model.geo.addPoint(x - r, y, 0, lc)
        p4 = gmsh.model.geo.addPoint(x, y - r, 0, lc)
        c1 = gmsh.model.geo.addCircleArc(p1, p_cen, p2)
        c2 = gmsh.model.geo.addCircleArc(p2, p_cen, p3)
        c3 = gmsh.model.geo.addCircleArc(p3, p_cen, p4)
        c4 = gmsh.model.geo.addCircleArc(p4, p_cen, p1)
        return gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
    end

    # --- 2. 孔洞布局 (自动对齐5腿结构) ---
    hole_loops = Int32[]
    step = W_leg + W_arch
    
    # A. 上排3个大孔
    push!(hole_loops, addHole(0.0, Y_hole_top, R_hole_large, lc))
    push!(hole_loops, addHole(-step, Y_hole_top, R_hole_large, lc))
    push!(hole_loops, addHole( step, Y_hole_top, R_hole_large, lc))

    # B. 下排4个小孔
    arch_offset = W_leg/2.0 + W_arch/2.0
    push!(hole_loops, addHole(-arch_offset, Y_hole_bot, R_hole_small, lc))
    push!(hole_loops, addHole( arch_offset, Y_hole_bot, R_hole_small, lc))
    push!(hole_loops, addHole(-(arch_offset + step), Y_hole_bot, R_hole_small, lc))
    push!(hole_loops, addHole( (arch_offset + step), Y_hole_bot, R_hole_small, lc))

    # --- 3. 外轮廓绘制 ---
    lines = Int32[] 
    
    x_leg3_right = W_leg / 2.0
    x_arch3_right = x_leg3_right + W_arch
    x_leg4_right = x_arch3_right + W_leg
    x_arch4_right = x_leg4_right + W_arch
    x_leg5_right = x_arch4_right + W_leg+5
    
    x_leg5_left  = -x_leg5_right
    x_arch1_left = -x_arch4_right
    x_leg4_left  = -x_leg4_right
    x_arch2_left = -x_arch3_right
    x_leg3_left  = -x_leg3_right
    x_arch3_right = -x_arch2_left # Type correction for coordinate consistency
    
    # --- A. 左侧飞翼 ---
    p_top_left = gmsh.model.geo.addPoint(-W_top_half, H_total, 0, lc)
    p_wing_left_tip = gmsh.model.geo.addPoint(x_leg5_left+8 , H_wing_tip, 0, lc)
    p_leg_start_top = gmsh.model.geo.addPoint(x_leg5_left+5, H_wing_start, 0, lc)
    p_leg_start_bot = gmsh.model.geo.addPoint(x_leg5_left+7, 0, 0, lc)

    Γ6 = gmsh.model.geo.addLine(p_top_left, p_wing_left_tip)
    Γ7 = gmsh.model.geo.addSpline([p_wing_left_tip,  p_leg_start_top])
    Γ8 = gmsh.model.geo.addLine(p_leg_start_top, p_leg_start_bot)
    push!(lines, Γ6, Γ7, Γ8)
    
    current_pt = p_leg_start_bot

    # --- B. 底部结构 ---
    # Leg 1
    p_l1_end = gmsh.model.geo.addPoint(x_arch1_left, 0, 0, lc)
    Γ1 = gmsh.model.geo.addLine(current_pt, p_l1_end)
    push!(lines, Γ1)
    
    # Arch 1
    p_a1_tl = gmsh.model.geo.addPoint(x_arch1_left, H_arch_straight, 0, lc)
    p_a1_cen = gmsh.model.geo.addPoint((x_arch1_left + x_leg4_left)/2, H_arch_straight, 0, lc)
    p_a1_tr = gmsh.model.geo.addPoint(x_leg4_left, H_arch_straight, 0, lc)
    p_a1_br = gmsh.model.geo.addPoint(x_leg4_left, 0, 0, lc)
    push!(lines, gmsh.model.geo.addLine(p_l1_end, p_a1_tl))
    push!(lines, gmsh.model.geo.addCircleArc(p_a1_tl, p_a1_cen, p_a1_tr, -1, 0, 0, -1)) 
    push!(lines, gmsh.model.geo.addLine(p_a1_tr, p_a1_br))

    # Leg 2
    p_l2_end = gmsh.model.geo.addPoint(x_arch2_left, 0, 0, lc)
    Γ2 = gmsh.model.geo.addLine(p_a1_br, p_l2_end)
    push!(lines, Γ2)

    # Arch 2
    p_a2_tl = gmsh.model.geo.addPoint(x_arch2_left, H_arch_straight, 0, lc)
    p_a2_cen = gmsh.model.geo.addPoint((x_arch2_left + x_leg3_left)/2, H_arch_straight, 0, lc)
    p_a2_tr = gmsh.model.geo.addPoint(x_leg3_left, H_arch_straight, 0, lc)
    p_a2_br = gmsh.model.geo.addPoint(x_leg3_left, 0, 0, lc)
    push!(lines, gmsh.model.geo.addLine(p_l2_end, p_a2_tl))
    push!(lines, gmsh.model.geo.addCircleArc(p_a2_tl, p_a2_cen, p_a2_tr, -1, 0, 0, -1))
    push!(lines, gmsh.model.geo.addLine(p_a2_tr, p_a2_br))

    # Leg 3
    p_l3_end = gmsh.model.geo.addPoint(x_leg3_right, 0, 0, lc)
    Γ3 = gmsh.model.geo.addLine(p_a2_br, p_l3_end)
    push!(lines, Γ3)

    # Arch 3
    p_a3_tl = gmsh.model.geo.addPoint(x_leg3_right, H_arch_straight, 0, lc)
    p_a3_cen = gmsh.model.geo.addPoint((x_leg3_right + x_arch3_right)/2, H_arch_straight, 0, lc)
    p_a3_tr = gmsh.model.geo.addPoint(x_arch3_right, H_arch_straight, 0, lc)
    p_a3_br = gmsh.model.geo.addPoint(x_arch3_right, 0, 0, lc)
    push!(lines, gmsh.model.geo.addLine(p_l3_end, p_a3_tl))
    push!(lines, gmsh.model.geo.addCircleArc(p_a3_tl, p_a3_cen, p_a3_tr, -1, 0, 0, -1))
    push!(lines, gmsh.model.geo.addLine(p_a3_tr, p_a3_br))

    # Leg 4
    p_l4_end = gmsh.model.geo.addPoint(x_leg4_right, 0, 0, lc)
    Γ4 = gmsh.model.geo.addLine(p_a3_br, p_l4_end)
    push!(lines, Γ4)

    # Arch 4
    p_a4_tl = gmsh.model.geo.addPoint(x_leg4_right, H_arch_straight, 0, lc)
    p_a4_cen = gmsh.model.geo.addPoint((x_leg4_right + x_arch4_right)/2, H_arch_straight, 0, lc)
    p_a4_tr = gmsh.model.geo.addPoint(x_arch4_right, H_arch_straight, 0, lc)
    p_a4_br = gmsh.model.geo.addPoint(x_arch4_right, 0, 0, lc)
    push!(lines, gmsh.model.geo.addLine(p_l4_end, p_a4_tl))
    push!(lines, gmsh.model.geo.addCircleArc(p_a4_tl, p_a4_cen, p_a4_tr, -1, 0, 0, -1))
    push!(lines, gmsh.model.geo.addLine(p_a4_tr, p_a4_br))

    # Leg 5
    p_l5_end = gmsh.model.geo.addPoint(x_leg5_right-7, 0, 0, lc)
    Γ5 = gmsh.model.geo.addLine(p_a4_br, p_l5_end)
    push!(lines, Γ5)
    
    # --- C. 右侧飞翼 ---
    p_leg_end_top = gmsh.model.geo.addPoint(x_leg5_right-5, H_wing_start, 0, lc)
    p_wing_right_tip = gmsh.model.geo.addPoint(x_leg5_right -8, H_wing_tip, 0, lc)
    p_top_right = gmsh.model.geo.addPoint(W_top_half, H_total, 0, lc)

    Γ9 = gmsh.model.geo.addLine(p_l5_end, p_leg_end_top)
    Γ10 = gmsh.model.geo.addSpline([p_leg_end_top,  p_wing_right_tip])
    Γ11 = gmsh.model.geo.addLine(p_wing_right_tip, p_top_right)
    push!(lines, Γ9, Γ10, Γ11)

    Γ12 = gmsh.model.geo.addLine(p_top_right, p_top_left)
    push!(lines, Γ12)

    # --- 4. 生成面 ---
    outer_loop = gmsh.model.geo.addCurveLoop(lines)
    surface_tags = [outer_loop; hole_loops]
    S1 = gmsh.model.geo.addPlaneSurface(surface_tags)
    gmsh.model.geo.synchronize()

    # 【修复重点 3】: 在 generateGeo 内部，必须使用 model.geo.addPhysicalGroup
    # 因为此时 geometry 还没 synchronize，model 中还不存在实体
    gmsh.model.geo.addPhysicalGroup(1, [Γ1], -1, "Γᵍ1")
    gmsh.model.geo.addPhysicalGroup(1, [Γ2], -1, "Γᵍ2")
    gmsh.model.geo.addPhysicalGroup(1, [Γ3], -1, "Γᵍ3")
    gmsh.model.geo.addPhysicalGroup(1, [Γ4], -1, "Γᵍ4")
    gmsh.model.geo.addPhysicalGroup(1, [Γ5], -1, "Γᵍ5")

    gmsh.model.geo.addPhysicalGroup(1, [Γ8], -1, "Γᵍl")
    gmsh.model.geo.addPhysicalGroup(1, [Γ9], -1, "Γᵍr")
    gmsh.model.geo.addPhysicalGroup(1, [Γ12], -1, "Γᵍt")
    gmsh.model.geo.addPhysicalGroup(1, [Γ6,Γ7], -1, "Γᵗl")
    gmsh.model.geo.addPhysicalGroup(1, [Γ11,Γ10], -1, "Γᵗr")
    gmsh.model.addPhysicalGroup(2, [S1], -1, "Ω")

    # 返回面 Tag 和 边界线列表 (以备后用)
    return S1, surface_tags
end

end