# Legacy reference implementations

- `lmp_swapstat4_main_bettercomments.cpp` is the original 4,590-line analysis
  program. It depends on in-house C++ headers that are not included here and is
  retained as a behavioral reference only.
- `lmp_swapstat4_main.py` is the former single-file core port. It is superseded
  by the package in `src/crumpling`.

Neither file is used by the current CLI or library API.
