# Legacy reference implementations

- `lmp_swapstat4_main_bettercomments.cpp` is the original 4,590-line analysis
  program. It depends on an in-house C++ math/stats library that is not included
  here.
- `lmp_swapstat4_main.py` is the former single-file core port of the detection and
  primary catalog writers.
- `cpp_core/` builds a **self-contained core binary** from the original `.cpp`:
  - `vendor/` — minimal stubs for `INIReader`, 3D vectors, `stoi`/`stod`, etc.
  - `build_core.py` — extracts parse/detect/catalog/structure sections (skips
    histogram statistics) and compiles `build/lmp_swapstat4_core`.

Neither the monolithic Python file nor the C++ core binary is used by the current
CLI; they are for parity checks (see `notebooks/compare_cpp_vs_python.ipynb`).
