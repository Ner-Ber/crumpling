#!/usr/bin/env bash
# Re-register the Jupyter kernel for the crumpling conda env.
# Run this if Cursor/Jupyter stops finding "Python (crumpling)".
set -euo pipefail

PYTHON="${CRUMPLING_PYTHON:-/home/neriberman/miniconda3/envs/crumpling/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
  echo "crumpling python not found: $PYTHON" >&2
  exit 1
fi

"$PYTHON" -m ipykernel install --user --name crumpling --display-name "Python (crumpling)"

KERNEL_JSON="$HOME/.local/share/jupyter/kernels/crumpling/kernel.json"
python3 - "$PYTHON" "$KERNEL_JSON" <<'PY'
import json, sys
from pathlib import Path
python, path = sys.argv[1], Path(sys.argv[2])
spec = {
    "argv": [python, "-Xfrozen_modules=off", "-m", "ipykernel_launcher", "-f", "{connection_file}"],
    "display_name": "Python (crumpling)",
    "language": "python",
    "metadata": {"debugger": True},
}
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(spec, indent=1) + "\n")
print(f"wrote {path}")
print(f"python={python}")
PY

"$PYTHON" -c "import crumpling, ipykernel; print('imports ok', crumpling.__file__)"
echo "Select notebook kernel: Python (crumpling)"
