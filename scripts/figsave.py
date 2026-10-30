#!/usr/bin/env python3
"""
figsave.py — one place where every figure is written, in every format the
project consumes. Import and call `save_fig(fig, "figures/figN_name.png")`
instead of `fig.savefig(...)`.

WHY THREE FORMATS
  .png  figures/           300 dpi raster. Required by python-pptx (the deck
                           builder cannot embed SVG or PDF).
  .pdf  figures/pdf/       TRUE VECTOR. Text stays text, lines stay lines —
                           scales without pixelation, and is what a journal
                           expects. Previously these were produced by
                           re-wrapping the PNG bitmap in a PDF container
                           (figs_to_pdf.py), which cannot recover vector
                           information discarded at rasterization.
  .svg  figures/svg/       TRUE VECTOR, for the HTML white paper. WeasyPrint
                           renders SVG natively but will NOT embed a PDF, so
                           SVG (not PDF) is the format that makes the white
                           paper's figures vector.

FONT EMBEDDING
  pdf.fonttype / ps.fonttype = 42 embeds TrueType fonts, so glyphs remain
  selectable, searchable and re-encodable. The matplotlib default (Type 3)
  is rejected by many journal production systems. Set explicitly here so the
  output does not depend on the installed matplotlib version.

CAVEAT — what vector output cannot fix
  Vector export preserves geometry that was drawn as geometry. Anything drawn
  as pixel data (imshow / pcolormesh / a rasterized inset) stays a bitmap
  inside the PDF; only its surrounding axes, ticks and labels become vector.
  All figures in this project are line/marker/patch art, so they vectorize
  fully — `verify_vector()` below is the check, not the assumption.
"""
import os
import matplotlib

# Must be set before any figure is written. Type 42 = embedded TrueType.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"   # keep text as <text>, not paths

DPI = 300
PDF_DIR = "figures/pdf"
SVG_DIR = "figures/svg"


def save_fig(fig, png_path, dpi=DPI, formats=("png", "pdf", "svg"), **kw):
    """Write `fig` to figures/<name>.png plus vector twins. Returns paths written.

    png_path is the canonical name (e.g. "figures/fig3_forest.png"); the vector
    outputs are derived from it, so a figure is named once.
    """
    kw.setdefault("bbox_inches", "tight")
    name = os.path.splitext(os.path.basename(png_path))[0]
    # Vector twins live in pdf/ and svg/ BESIDE the png, derived from its own folder.
    # For the usual "figures/x.png" this is figures/pdf and figures/svg exactly as
    # before; for a sub-collection like "figures/bolte/x.png" it becomes
    # figures/bolte/{pdf,svg}, so a parallel analysis keeps its outputs together
    # instead of mixing its vectors into the main figure set.
    base = os.path.dirname(png_path) or "."
    pdf_dir = PDF_DIR if base == "figures" else os.path.join(base, "pdf")
    svg_dir = SVG_DIR if base == "figures" else os.path.join(base, "svg")
    written = []
    if "png" in formats:
        os.makedirs(base, exist_ok=True)
        fig.savefig(png_path, dpi=dpi, **kw)
        written.append(png_path)
    if "pdf" in formats:
        os.makedirs(pdf_dir, exist_ok=True)
        p = os.path.join(pdf_dir, name + ".pdf")
        fig.savefig(p, **kw)            # no dpi -> vector, not rasterized
        written.append(p)
    if "svg" in formats:
        os.makedirs(svg_dir, exist_ok=True)
        p = os.path.join(svg_dir, name + ".svg")
        fig.savefig(p, **kw)
        written.append(p)
    return written


def verify_vector(pdf_path):
    """Inspect a PDF's content streams: is the drawing vector or a pasted bitmap?

    Returns {images, fonts, vector_ops, text_ops, is_vector}. A PDF produced by
    wrapping a PNG has images>=1 and vector_ops==0; a true matplotlib vector PDF
    has thousands of path operators and real text-showing operators.
    """
    import re, zlib
    b = open(pdf_path, "rb").read()
    n_img = len(re.findall(rb"/Subtype\s*/Image", b))
    n_font = len(re.findall(rb"/Type\s*/Font", b))
    ops = txt = 0
    for m in re.finditer(rb"stream\r?\n", b):
        s = b[m.end():]
        e = s.find(b"endstream")
        if e < 0:
            continue
        try:
            d = zlib.decompress(s[:e])
        except Exception:
            continue
        ops += len(re.findall(rb"\b(?:re|m|l|c)\b", d))
        txt += len(re.findall(rb"\b(?:Tj|TJ)\b", d))
    return {"images": n_img, "fonts": n_font, "vector_ops": ops,
            "text_ops": txt, "is_vector": ops > 0 and txt > 0}
