#!/usr/bin/env python3
"""Place a publication-resolution PNG on an exactly sized PDF page."""

from pathlib import Path
import sys

from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(
            "usage: raster_png_to_pdf.py INPUT.png OUTPUT.pdf WIDTH_IN HEIGHT_IN"
        )

    input_path = Path(sys.argv[1]).resolve()
    output_path = Path(sys.argv[2]).resolve()
    width_pt = float(sys.argv[3]) * 72.0
    height_pt = float(sys.argv[4]) * 72.0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(
        str(output_path),
        pagesize=(width_pt, height_pt),
        pageCompression=1,
    )
    pdf.drawImage(
        ImageReader(str(input_path)),
        0,
        0,
        width=width_pt,
        height=height_pt,
        preserveAspectRatio=False,
        mask="auto",
    )
    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    main()
