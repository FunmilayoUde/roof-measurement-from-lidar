"""Build a standalone copy of the write-up for any static host (GitHub Pages, Netlify, ...).

    uv run python scripts/build_site.py

docs/report/index.html is a page fragment without <html>/<head> tags;
this wraps it in a full document and copies report_data.json beside it, into site/.
"""

from __future__ import annotations

import shutil
from pathlib import Path

SRC = Path("docs/report")
OUT = Path("site")

HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="description" content="Roof surfaces, pitch, true area and typed edges measured from aerial LiDAR alone
(Dutch AHN), with an interactive 3D model.">
<style>body { margin: 0; } img { max-width: 100%; } [hidden] { display: none !important; }</style>
"""


def main() -> None:
    page = (SRC / "index.html").read_text()
    head_end = page.index("</style>") + len("</style>")
    OUT.mkdir(exist_ok=True)
    (OUT / "index.html").write_text(
        HEAD + page[:head_end] + "\n</head>\n<body>\n" + page[head_end:] + "\n</body>\n</html>\n"
    )
    shutil.copy(SRC / "report_data.json", OUT / "report_data.json")
    print(f"wrote {OUT}/index.html and {OUT}/report_data.json")


if __name__ == "__main__":
    main()
