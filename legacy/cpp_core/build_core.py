#!/usr/bin/env python3
"""Build a self-contained core-only binary from the legacy C++ reference.

The original program depends on an in-house math/stats library that is not
shipped in this repository. This script:

1. Extracts the core pipeline sections (parse → detect → write catalogs /
   avalanche structure), skipping histogram/statistics blocks.
2. Rewrites includes to the local stubs under ``vendor/``.
3. Compiles ``build/lmp_swapstat4_core`` with g++.

Usage:
    python build_core.py
"""

from __future__ import annotations

import pathlib
import subprocess
import sys


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SRC = ROOT / "legacy" / "lmp_swapstat4_main_bettercomments.cpp"
BUILD = HERE / "build"
OUT_CPP = BUILD / "lmp_swapstat4_main_core.cpp"
OUT_BIN = BUILD / "lmp_swapstat4_core"
VENDOR = HERE / "vendor"

# 1-based inclusive line ranges kept from the original source.
KEEP_RANGES = (
    (1, 1746),    # through Detect avalanches
    (3257, 3402), # Write markup / subav / av / all-av catalogs
    (4351, 4584), # Write avalanche structure
)


HEADER = """\
// Auto-generated core-only build of lmp_swapstat4_main_bettercomments.cpp
// Histogram/statistics sections are omitted; catalogs + structure only.
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <algorithm>
#include <map>
#include <cmath>
#include <sstream>

#include "general/stringext.h"
#include "aside/inih/cpp/INIReader.h"
#include "math/basic_math/basic_math.h"
#include "math/geometry/planar/vector.h"
#include "math/geometry/planar/rect.h"
#include "math/geometry/planar/circle.h"
#include "phys/atomic/pathanalysis.h"

using namespace std;

"""


def extract_core(source_text: str) -> str:
    lines = source_text.splitlines(keepends=True)
    chunks: list[str] = [HEADER]
    # Absolute original line numbers; first range drops the private includes / using.
    for index, (start, end) in enumerate(KEEP_RANGES):
        slice_start = 29 if start == 1 else start  # structs begin after `using namespace std;`
        chunks.append(f"\n/* ---- original lines {slice_start}-{end} ---- */\n")
        chunks.extend(lines[slice_start - 1 : end])
        if index == 0:
            chunks.append(
                "\n"
                "    cout << \"\\n[core build] skipping histogram/statistics sections...\\n\";\n"
            )
    chunks.append("\n    return 0;\n}\n")
    return "".join(chunks)


def compile_binary(cpp_path: pathlib.Path, bin_path: pathlib.Path) -> None:
    cmd = [
        "g++",
        "-std=c++17",
        "-O2",
        "-Wall",
        "-Wno-unused-variable",
        "-Wno-unused-but-set-variable",
        "-Wno-sign-compare",
        f"-I{VENDOR}",
        str(cpp_path),
        "-o",
        str(bin_path),
        "-lm",
    ]
    print("compile:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> int:
    if not SRC.is_file():
        print(f"missing source: {SRC}", file=sys.stderr)
        return 1
    BUILD.mkdir(parents=True, exist_ok=True)
    core = extract_core(SRC.read_text(encoding="utf-8", errors="replace"))
    OUT_CPP.write_text(core, encoding="utf-8")
    print(f"wrote {OUT_CPP} ({len(core.splitlines())} lines)")
    compile_binary(OUT_CPP, OUT_BIN)
    print(f"binary: {OUT_BIN}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
