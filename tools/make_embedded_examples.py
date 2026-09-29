"""Writes Examples/embedded/*.js: gzip + base64 copies of the example files.

A page opened from disk (file://) may not fetch Examples/*, but it can load
scripts next to it, so the example buttons read these copies instead. Each
script sets window.__3dsExamples["Examples/<path>"] = "<base64 of gzip>".
PMN_300k.rmc6f (30 MB, 13 MB as base64) is left out; from disk its button
opens the file dialog.

Usage (from the repository root, after changing an example):
    python tools/make_embedded_examples.py
"""
import base64
import gzip
from pathlib import Path

FILES = [
    "Examples/LiFeO2.rmc6f",
    "Examples/CaTiO3.rmc6f",
    "Examples/magnetic/MnO_order_spins.txt",
    "Examples/magnetic/MnO_random_spins.txt",
    "Examples/magnetic/MnO.rmc6f",
    "Examples/magnetic/MnO_spins.cfg",
    "Examples/magnetic/MnO.dat",
]
OUT = Path("Examples/embedded")


def embedded_name(path):
    return path[len("Examples/"):].replace("/", "_") + ".js"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for path in FILES:
        data = Path(path).read_bytes()
        # mtime=0 keeps the output identical from run to run
        packed = base64.b64encode(gzip.compress(data, compresslevel=9, mtime=0)).decode("ascii")
        script = f'(window.__3dsExamples = window.__3dsExamples || {{}})["{path}"] = "{packed}";\n'
        (OUT / embedded_name(path)).write_text(script, encoding="ascii", newline="\n")
        print(f"{path}: {len(data)} -> {len(packed)} bytes")


if __name__ == "__main__":
    main()
