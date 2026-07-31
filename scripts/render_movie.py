#!/usr/bin/env python3
# ============================================================================
#  render_movie.py -- Animation 3D de la deformation de la bulle
# ----------------------------------------------------------------------------
#  Strategie HPC : la DNS ne fait pas de graphique (cout), elle ecrit des
#  snapshots "snapshot-XXXXXX.XX". Ce script rend chaque snapshot en PNG via
#  l'outil Basilisk ./render (make render), puis assemble une video avec ffmpeg.
#
#  Prerequis :
#    make render          # compile bin ./render (OSMesa/GL requis)
#    ffmpeg dans le PATH
#
#  Usage :
#    python scripts/render_movie.py --dir simulations/lvl9
#    python scripts/render_movie.py --dir simulations/lvl9 --fps 20 --out bulle.mp4
# ============================================================================
import argparse
import glob
import os
import re
import subprocess
import sys


def natural_key(path):
    m = re.search(r"snapshot-([0-9.]+)", os.path.basename(path))
    return float(m.group(1)) if m else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="dossier de simulation")
    ap.add_argument("--render", default="./render", help="binaire de rendu Basilisk")
    ap.add_argument("--out", default="bubble.mp4")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--pattern", default="snapshot-*")
    args = ap.parse_args()

    sim = os.path.abspath(args.dir)
    render = os.path.abspath(args.render)
    if not os.path.isfile(render):
        sys.exit(f"Binaire de rendu absent : {render}\n  -> lancez 'make render'")

    snaps = sorted(glob.glob(os.path.join(sim, args.pattern)), key=natural_key)
    snaps = [s for s in snaps if not s.endswith(".png")]
    if not snaps:
        sys.exit(f"Aucun snapshot ({args.pattern}) dans {sim}")

    frames_dir = os.path.join(sim, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    print(f"  {len(snaps)} snapshots a rendre...")

    for k, snap in enumerate(snaps):
        png = os.path.join(frames_dir, f"frame_{k:04d}.png")
        print(f"  [{k+1}/{len(snaps)}] {os.path.basename(snap)} -> {os.path.basename(png)}")
        r = subprocess.run([render, snap, png],
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if r.returncode != 0:
            print("    erreur de rendu :", r.stderr.decode(errors="ignore"))

    out = os.path.join(sim, args.out)
    cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
           "-i", os.path.join(frames_dir, "frame_%04d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", out]
    print("  assemblage :", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"  Video : {out}")


if __name__ == "__main__":
    main()
