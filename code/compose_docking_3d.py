from PIL import Image, ImageDraw, ImageFont
import os

BASE = r"D:\workspace\research\protein-drug\06-figures\main\pymol"
OUT = r"D:\workspace\research\protein-drug\06-figures\main"
whole = Image.open(os.path.join(BASE, "whole.png")).convert("RGB")
pocket = Image.open(os.path.join(BASE, "pocket_clean.png")).convert("RGB")

H = 1400
def fit(im, h=H):
    w = int(im.width * h / im.height)
    return im.resize((w, h), Image.LANCZOS)

whole, pocket = fit(whole), fit(pocket)
pad = 30
W = whole.width + pocket.width + pad * 3
canvas = Image.new("RGB", (W, H + pad * 2), "white")
d = ImageDraw.Draw(canvas)
font = ImageFont.truetype(r"C:\Windows\Fonts\arialbd.ttf", 64)

canvas.paste(whole, (pad, pad))
canvas.paste(pocket, (pad * 2 + whole.width, pad))
d.text((pad + 18, pad + 10), "a", fill="black", font=font)
d.text((pad * 2 + whole.width + 18, pad + 10), "b", fill="black", font=font)

out_png = os.path.join(OUT, "fig_ga_docking_3d.png")
out_pdf = os.path.join(OUT, "fig_ga_docking_3d.pdf")
canvas.save(out_png)
canvas.save(out_pdf, "PDF", resolution=300)
print("saved", out_png, canvas.size)
