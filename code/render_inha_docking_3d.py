# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#   "pymol-open-source-whl",
# ]
# ///
"""Render the entrectinib--4DRE complex and its binding pocket with PyMOL/OSMesa."""
import os

os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol  # pytype: disable=import-error

pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd  # pytype: disable=import-error

PDB = r"D:\scratch\protein-drug\md\entrectinib_4DRE\complex_H.pdb"
OUT = r"D:\workspace\research\protein-drug\06-figures\main\pymol"
os.makedirs(OUT, exist_ok=True)

cmd.load(PDB, "complex")
n = cmd.count_atoms("all")
print("atoms loaded:", n)
if n == 0:
    print("ERROR: no atoms loaded")
    cmd.quit()
    raise SystemExit(1)

cmd.remove("hydro")
cmd.hide("everything")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("cartoon_transparency", 0.10)
cmd.set("cartoon_side_chain_helper", 1)
cmd.set("stick_radius", 0.16)
cmd.set("surface_quality", 1)
cmd.set("depth_cue", 0)

# protein
cmd.show("cartoon", "polymer")
cmd.color("gray75", "polymer")
cmd.set("cartoon_color", "gray75", "polymer")

# cofactor and ligand
cmd.show("sticks", "resn NAI")
cmd.color("teal", "resn NAI and elem C")
cmd.color("red", "resn NAI and elem O")
cmd.color("blue", "resn NAI and elem N")
cmd.color("orange", "resn NAI and elem P")

cmd.show("sticks", "resn LIG")
cmd.color("orange", "resn LIG and elem C")
cmd.color("red", "resn LIG and elem O")
cmd.color("blue", "resn LIG and elem N")
cmd.color("gray30", "resn LIG and elem F")

# key pocket residues
key = "resi 16+96+98+149+158 and polymer"
cmd.select("keyres", key)
cmd.show("sticks", "keyres")
cmd.color("green", "keyres and elem C")
cmd.color("red", "keyres and elem O")
cmd.color("blue", "keyres and elem N")
cmd.set("cartoon_transparency", 0.25, "keyres")

# whole complex
cmd.orient("polymer")
cmd.zoom("polymer", 2.0)
cmd.png(os.path.join(OUT, "whole.png"), width=1600, height=1400, dpi=300)

# pocket view
cmd.hide("cartoon", "polymer")
cmd.show("cartoon", "polymer")
cmd.show("surface", "polymer")
cmd.set("transparency", 0.78)
cmd.zoom("resn LIG", 7.0)
cmd.turn("y", 25)
cmd.turn("x", 12)
cmd.png(os.path.join(OUT, "pocket_clean.png"), width=1600, height=1400, dpi=300)

cmd.save(os.path.join(OUT, "inha_entrectinib.pse"))
print("saved whole_clean.png, pocket_clean.png, inha_entrectinib.pse")
cmd.quit()
