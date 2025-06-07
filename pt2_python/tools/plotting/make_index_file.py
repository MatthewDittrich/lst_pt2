import os
import warnings
from pathlib import Path

def generate_html_indexes(base_dir):
    base_dir = Path(base_dir)

    for root, dirs, files in os.walk(base_dir):
        root_path = Path(root)

        # Find matching .png and .pdf pairs
        pngs = {f.stem for f in root_path.glob("*.png")}
        pdfs = {f.stem for f in root_path.glob("*.pdf")}
        common = sorted(pngs & pdfs)

        if not common:
            continue  # Skip if no valid pairs

        # Warn about unmatched files
        unmatched = (pngs ^ pdfs)
        for name in unmatched:
            warnings.warn(f"Missing pair for: {name} in {root_path}")

        # Start writing index.html
        index_path = root_path / "index.html"
        with open(index_path, "w") as f:
            f.write("<!DOCTYPE html>\n<html>\n<head>\n")
            f.write("<meta charset='utf-8'>\n")
            f.write("<title>Plot Index</title>\n")
            f.write("""
<style>
body { font-family: sans-serif; margin: 20px; }
.grid { display: flex; flex-wrap: wrap; gap: 20px; }
.item { width: calc(25% - 20px); box-shadow: 0 0 5px #ccc; padding: 10px; background: #f9f9f9; }
.item img { width: 100%; height: auto; border: 1px solid #ccc; }
.subdirs { margin-bottom: 30px; }
.subdirs a { display: inline-block; margin-right: 15px; font-weight: bold; }
</style>
</head>
<body>
""")
            f.write(f"<h1>Index of {root_path.relative_to(base_dir)}</h1>\n")

            # Subdirectory links
            f.write("<div class='subdirs'>\n")
            for subdir in sorted(dirs):
                sub_index = root_path / subdir / "index.html"
                if sub_index.exists():
                    f.write(f"<a href='{subdir}/index.html'>{subdir}</a>\n")
            f.write("</div>\n")

            # Plot grid
            f.write("<div class='grid'>\n")
            for name in common:
                png_file = f"{name}.png"
                pdf_file = f"{name}.pdf"
                f.write(f"<div class='item'>\n")
                f.write(f"<a href='{pdf_file}'><img src='{png_file}' alt='{name}'></a>\n")
                f.write(f"<p>{name}</p>\n")
                f.write("</div>\n")
            f.write("</div>\n</body>\n</html>")

        print(f"Created {index_path}")
