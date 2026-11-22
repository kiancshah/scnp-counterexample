"""Figures: competing chains of a bad interval, drawn by rank, as a spiral, and as a rosette.

Usage:  python3 figures.py OUTDIR
For a Boolean interval whose competing chains are rotations of one cyclic order a_1..a_k of the atoms, an element v
is placed at the complex number  z(v) = sum of exp(2 pi i j / k) over the atoms a_j below v,  at height rank(v).
"""
import sys, math, cmath
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from core import *
from spiral import structure, chains_with_g, atom_orders

COL = ["#3b6fd4", "#d1582c", "#1c9a78", "#c7457b", "#a97a12", "#7a55c7", "#4f8d1a", "#c23a3a", "#2a8c9c", "#555555"]
GREY = "#c4cccb"; INK = "#15211f"; MUTE = "#7d8b88"; RING = "#b23a2e"
plt.rcParams["font.family"] = "DejaVu Sans"

def name(sw_):
    return "(" + ",".join(str(x) for x in sw_) + ")"

def apply_steps(u, steps):
    out = [u]; v = list(u)
    for a, b in steps:
        v[a - 1], v[b - 1] = v[b - 1], v[a - 1]; out.append(tuple(v))
    return out

def cover_edges(st, w):
    E = []
    for lev in st["levels"][:-1]:
        for v in lev:
            for a, b, x in covers_up(v):
                if x in st["below"]:
                    E.append((v, x))
    return E

def hasse(us, ws, fn, width=7.6, height=5.4):
    u, w = P(us), P(ws)
    st = structure(u, w); mins = minimal_g(u, w)
    pos = {}
    L = st["L"]
    for r, lev in enumerate(st["levels"]):
        srt = sorted(lev, reverse=True); k = len(srt)
        sp = min(10.0 / (k - 1), 1.3) if k > 1 else 0
        for j, v in enumerate(srt):
            pos[v] = (-sp * (k - 1) / 2 + sp * j, r)
    fig = plt.figure(figsize=(width, height)); ax = fig.add_axes([0, 0, 1, 1])
    ax.add_collection(LineCollection([[pos[a], pos[b]] for a, b in cover_edges(st, w)], colors=GREY, linewidths=0.5, zorder=1))
    used = set()
    for i, m in enumerate(mins):
        chs = chains_with_g(u, w, m)
        segs = set()
        for ch in chs:
            seq = apply_steps(u, ch)
            for a, b in zip(seq, seq[1:]):
                segs.add((a, b)); used.add(a); used.add(b)
        off = (i - (len(mins) - 1) / 2) * 0.035
        for a, b in segs:
            ax.plot([pos[a][0] + off, pos[b][0] + off], [pos[a][1], pos[b][1]], color=COL[i], lw=1.5, solid_capstyle="round", zorder=3)
        first = sorted({ch[0] for ch in chs})
        ax.plot([], [], color=COL[i], lw=2, label="first step " + ", ".join(name(s) for s in first) + f"  ({len(chs)} chains)")
    xs = [pos[v][0] for v in st["nodes"]]; ys = [pos[v][1] for v in st["nodes"]]
    ax.scatter(xs, ys, s=7, color=[INK if v in used else "#aab4b2" for v in st["nodes"]], zorder=4, linewidths=0)
    ax.text(pos[u][0], pos[u][1] - 0.32, "u = " + S(u), ha="center", va="top", fontsize=9, color=INK)
    ax.text(pos[w][0], pos[w][1] + 0.28, "w = " + S(w), ha="center", va="bottom", fontsize=9, color=INK)
    ax.legend(loc="lower center", fontsize=9, frameon=False, ncol=3, columnspacing=1.4, handlelength=1.6)
    ax.set_xlim(-5.6, 5.6); ax.set_ylim(-2.1, L + 0.8); ax.axis("off")
    fig.savefig(fn, bbox_inches="tight", pad_inches=0.03); plt.close(fig)

def cube_data(us, ws):
    u, w = P(us), P(ws)
    st = structure(u, w); assert st["boolean"]
    mins = minimal_g(u, w)
    chains = [chains_with_g(u, w, m) for m in mins]
    orders = [atom_orders(u, w, st, c[0]) for c in chains]
    base = orders[0]; k = len(base)
    rot = {tuple(base[i:] + base[:i]): i for i in range(k)}
    for cs in chains:
        for c in cs:
            assert tuple(atom_orders(u, w, st, c)) in rot
    start = {rot[tuple(o)]: i for i, o in enumerate(orders)}
    inv_swap = {v: a for a, v in st["atom_swap"].items()}
    ang = {inv_swap[s]: 2 * math.pi * j / k for j, s in enumerate(base)}
    z = {v: sum(cmath.exp(1j * ang[a]) for a in st["below"][v]) for v in st["nodes"]}
    rank = {v: r for r, lev in enumerate(st["levels"]) for v in lev}
    return dict(u=u, w=w, st=st, base=base, k=k, start=start, z=z, rank=rank, ang=ang, inv_swap=inv_swap)

def rotation_points(D, r):
    k = D["k"]; pts = [0j]; zz = 0j
    for j in range(k):
        zz += cmath.exp(2j * math.pi * ((r + j) % k) / k); pts.append(zz)
    return pts

def spiral(us, ws, fn, yaw=0.55, pitch=0.30, width=5.2, height=5.6):
    D = cube_data(us, ws); st = D["st"]; k = D["k"]; L = st["L"]
    R = 1 / math.sin(math.pi / k)
    sx = 1.0; sh = 2.0 * R / L * 1.35
    cy, sy = math.cos(yaw), math.sin(yaw); cp, sp = math.cos(pitch), math.sin(pitch)
    def proj(zc, h):
        X = zc.real * cy - zc.imag * sy; Z = zc.real * sy + zc.imag * cy
        return (X * sx, h * sh * cp - Z * sp)
    pr = {v: proj(D["z"][v], D["rank"][v]) for v in st["nodes"]}
    fig = plt.figure(figsize=(width, height)); ax = fig.add_axes([0, 0, 1, 1])
    big = len(st["nodes"]) > 300
    ax.add_collection(LineCollection([[pr[a], pr[b]] for a, b in cover_edges(st, D["w"])], colors=GREY,
                                     linewidths=0.18 if big else 0.45, alpha=0.55 if big else 0.9, zorder=1))
    for r in range(k):
        if r not in D["start"]: continue
        pts = rotation_points(D, r); c = COL[D["start"][r] % len(COL)]
        xy = [proj(p, h) for h, p in enumerate(pts)]
        ax.plot([p[0] for p in xy], [p[1] for p in xy], color=c, lw=1.9, solid_joinstyle="round", zorder=3)
        ax.scatter([p[0] for p in xy[1:-1]], [p[1] for p in xy[1:-1]], s=9, color=c, zorder=4, linewidths=0)
    ring = [proj(cmath.exp(2j * math.pi * j / k), 1) for j in range(k + 1)]
    ax.plot([p[0] for p in ring], [p[1] for p in ring], color=RING, lw=1.3, zorder=5)
    bot = proj(0j, 0); top = proj(D["z"][D["w"]], L)
    ax.scatter([bot[0], top[0]], [bot[1], top[1]], s=22, color=INK, zorder=6)
    for j, s in enumerate(D["base"]):
        q = proj((1.0 + 0.5 * R) * cmath.exp(2j * math.pi * j / k), 1)
        ax.text(q[0], q[1], name(s), ha="center", va="center", fontsize=9, color=INK, zorder=7)
    ax.text(bot[0], bot[1] - 0.22 * R, "u", ha="center", va="top", fontsize=10, color=INK)
    ax.text(top[0], top[1] + 0.12 * R, "w", ha="center", va="bottom", fontsize=10, color=INK)
    xs = [p[0] for p in pr.values()]; ys = [p[1] for p in pr.values()]
    m = 0.45 * R
    ax.set_xlim(min(xs) - m, max(xs) + m); ax.set_ylim(min(ys) - 1.4 * m, max(ys) + m); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(fn, bbox_inches="tight", pad_inches=0.03); plt.close(fig)

def rosette(us, ws, fn, width=5.2, height=5.2):
    D = cube_data(us, ws); st = D["st"]; k = D["k"]
    R = 1 / math.sin(math.pi / k)
    fig = plt.figure(figsize=(width, height)); ax = fig.add_axes([0, 0, 1, 1])
    big = len(st["nodes"]) > 300
    segs = [[(D["z"][a].real, D["z"][a].imag), (D["z"][b].real, D["z"][b].imag)] for a, b in cover_edges(st, D["w"])]
    ax.add_collection(LineCollection(segs, colors=GREY, linewidths=0.18 if big else 0.4, alpha=0.5, zorder=1))
    for r in range(k):
        pts = rotation_points(D, r); xs = [p.real for p in pts]; ys = [p.imag for p in pts]
        mid = pts[k // 2]; d = abs(mid); lab = name(D["base"][r])
        lp = mid / d * (d + 0.20 * R + 0.40)
        if r in D["start"]:
            c = COL[D["start"][r] % len(COL)]
            ax.plot(xs, ys, color=c, lw=1.9, zorder=3, solid_joinstyle="round")
            ax.scatter(xs[1:-1], ys[1:-1], s=9, color=c, zorder=4, linewidths=0)
            ax.text(lp.real, lp.imag, lab, color=c, fontsize=10, ha="center", va="center", zorder=6)
        else:
            ax.plot(xs, ys, color=MUTE, lw=1.1, ls=(0, (3, 3)), zorder=2)
            ax.text(lp.real, lp.imag, lab + "\nmissing", color=MUTE, fontsize=9, ha="center", va="center", zorder=6)
    ring = [cmath.exp(2j * math.pi * j / k) for j in range(k + 1)]
    ax.plot([p.real for p in ring], [p.imag for p in ring], color=RING, lw=1.3, zorder=5)
    ax.scatter([0], [0], s=26, color=INK, zorder=7)
    lim = 2 * R + 0.16 * R + 0.75
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(fn, bbox_inches="tight", pad_inches=0.03); plt.close(fig)

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    hasse("214365", "541632", out + "/hasse_s6.pdf")
    spiral("1432765", "5276143", out + "/spiral_s7.pdf"); rosette("1432765", "5276143", out + "/rosette_s7.pdf")
    spiral("3217654A98", "5482A19376", out + "/spiral_s10.pdf"); rosette("3217654A98", "5482A19376", out + "/rosette_s10.pdf")
    spiral("1324", "4231", out + "/spiral_s4.pdf"); rosette("1324", "4231", out + "/rosette_s4.pdf")
    print("figures written to", out)
