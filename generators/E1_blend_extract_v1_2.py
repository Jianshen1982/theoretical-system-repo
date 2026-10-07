# -*- coding: utf-8 -*-
# ============================================================
# E1_blend_extract_v1_2.py  筛查线(Blender)—— 球心顶点直读导出
# 纪律: D.0 —— 本脚本输出不进入定案链(Blender 顶点底层 float32)
# 补正: δ空间/δ镜像空间 = 4 顶点四面体可视化辅助结构(非测量主体);
#       测量主体 = "+Z闭合构型"(13 顶点球心模型)
# 运行: blender --background "D:\物理数学\数学\理想球体的有限不可实现性\+Z闭合构型.blend" \
#            --python E1_blend_extract_v1_2.py
# ============================================================
import bpy, json, sys
from datetime import datetime
from pathlib import Path

OUT_ROOT  = Path(r"D:\理论体系\repo\text")
CODE, VER = "E1_blend_extract", "v1_2"
# role: primary=测量主体(13 球心) / aux=可视化辅助(4 顶点, 不入定案)
OBJECTS = {
    "+Z闭合构型":         dict(role="primary", expect=13),
    "+Z闭合构型δ空间":     dict(role="aux",     expect=4),
    "+Z闭合构型δ镜像空间":  dict(role="aux",     expect=4),
}

def fail(msg):
    print("[FAIL-FAST]", msg); sys.exit(1)

def out_path(code, ver):
    now = datetime.now()
    day = OUT_ROOT / now.strftime("%Y") / now.strftime("%m") / now.strftime("%d")
    day.mkdir(parents=True, exist_ok=True)
    seq = 1 + len(list(day.glob(f"{code}_{ver}_*")))
    return day / f"{code}_{ver}_{seq:03d}_{now.strftime('%Y%m%d-%H%M%S')}.json"

def main():
    records = []
    for name, spec in OBJECTS.items():
        obj = bpy.data.objects.get(name)
        if obj is None: fail(f"缺少对象: {name}")
        M = obj.matrix_world
        pts = [tuple(float(c) for c in (M @ v.co)) for v in obj.data.vertices]
        flags = []
        if len(pts) != spec["expect"]:
            flags.append(f"顶点数 {len(pts)} ≠ 预期 {spec['expect']}"
                         + ("(球心顶点模型)" if spec["role"]=="primary" else "(四面体辅助结构)"))
        records.append(dict(object=name, role=spec["role"],
                            n_vertices=len(pts), vertices=pts, flags=flags))
        print(f"[E1] {name}[{spec['role']}]: 顶点 {len(pts)}  flags={flags}")
    payload = dict(code=CODE, ver=VER,
                   model=r"D:\物理数学\数学\理想球体的有限不可实现性\+Z闭合构型.blend",
                   blender=bpy.app.version_string,
                   precision_note="Blender 顶点底层 float32——筛查线数据,不入定案链(D.0)",
                   model_type="球心顶点模型;δ空间/δ镜像空间为4顶点四面体可视化辅助结构(非测量主体)",
                   created=datetime.now().isoformat(timespec="seconds"),
                   records=records)
    p = out_path(CODE, VER)
    p.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[E1] 导出:", p)

main()
