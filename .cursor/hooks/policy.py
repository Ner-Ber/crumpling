#!/usr/bin/env python3
"""Workspace hooks for change checks and new-module approval.

Limits (keep the rules in ``.cursor/rules`` in sync):

- ``REPAIR_LIMIT``: 3 failed checks, then stop and report
- ``BUDGET_SECONDS``: 10 minutes from the first failed check
- ``COMMAND_TIMEOUT_SECONDS``: 3 minutes for one automatic check
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPAIR_LIMIT = 3
BUDGET_SECONDS = 10 * 60
COMMAND_TIMEOUT_SECONDS = 170
MODERATE_CHARS = 500
MODERATE_FILES = 2
FOLLOWUP_MARK = "[verify-hook]"
OUTPUT_LIMIT = 4000

MODULE_TESTS = (
    ("crumpling/parsing.py", "tests/test_parsing.py"),
    ("crumpling/models.py", "tests/test_parsing.py"),
    ("crumpling/geometry.py", "tests/test_cascades.py"),
    ("crumpling/detection/cascades.py", "tests/test_cascades.py"),
    ("crumpling/detection/avalanches.py", "tests/test_avalanches.py"),
    ("crumpling/detection/", "tests/test_cascades.py"),
    ("crumpling/pipeline.py", "tests/test_pipeline.py"),
    ("crumpling/cli.py", "tests/test_pipeline.py"),
    ("crumpling/config.py", "tests/test_pipeline.py"),
    ("crumpling/__main__.py", "tests/test_pipeline.py"),
    ("crumpling/catalogs.py", "tests/test_catalogs.py"),
    ("crumpling/output/", "tests/test_catalogs.py"),
)

EXISTING_MODULES = (
    "crumpling/parsing.py",
    "crumpling/models.py",
    "crumpling/geometry.py",
    "crumpling/config.py",
    "crumpling/detection/cascades.py",
    "crumpling/detection/avalanches.py",
    "crumpling/pipeline.py",
    "crumpling/cli.py",
    "crumpling/catalogs.py",
    "crumpling/output/catalogs.py",
    "crumpling/output/structure.py",
)

APPROVE_RE = re.compile(
    r"\b(yes|yeah|yep|approve|approved|go ahead|do it|create it|ok|okay|lgtm)\b",
    re.IGNORECASE,
)
REJECT_RE = re.compile(
    r"\b(no|don't|do not|dont|wait|stop|hold on|not yet|never)\b",
    re.IGNORECASE,
)
PY_PATH_RE = re.compile(r"crumpling(?:/[\w-]+)+\.py")
DIR_PATH_RE = re.compile(r"crumpling(?:/[\w-]+)+")
CREATE_VERB_RE = re.compile(
    r"\b(touch|tee|mkdir|cp|mv|install|ln|open)\b|>",
)


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    raw = sys.stdin.read()
    try:
        payload = json.loads(raw) if raw.strip() else {}
    except json.JSONDecodeError:
        payload = {}
    if not isinstance(payload, dict):
        payload = {}
    try:
        result = dispatch(mode, payload)
    except Exception as exc:
        sys.stderr.write("policy hook failed: {0}\n".format(exc))
        result = {"permission": "allow"} if mode.startswith("guard") else {}
        if mode == "user-prompt":
            result = {"continue": True}
    sys.stdout.write(json.dumps(result))


def dispatch(mode: str, payload: Mapping[str, Any]) -> Dict[str, Any]:
    if mode == "edit":
        return handle_edit(payload)
    if mode == "shell-done":
        return handle_shell(payload, failed=False)
    if mode == "shell-fail":
        return handle_shell(payload, failed=True)
    if mode == "stop":
        return handle_stop(payload)
    if mode == "guard-write":
        return handle_guard_write(payload)
    if mode == "guard-shell":
        return handle_guard_shell(payload)
    if mode == "user-prompt":
        return handle_user_prompt(payload)
    return {}


def workspace_root(payload: Mapping[str, Any]) -> Path:
    roots = payload.get("workspace_roots") or []
    if isinstance(roots, list) and roots:
        return Path(str(roots[0])).resolve()
    return Path(__file__).resolve().parents[2]


def state_dir(root: Path) -> Path:
    return root / ".cursor" / "hooks" / "state"


def state_path(root: Path, conversation_id: str) -> Path:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", conversation_id)[:80] or "unknown"
    return state_dir(root) / (cleaned + ".json")


def default_state() -> Dict[str, Any]:
    return {
        "files": {},
        "verified": False,
        "repairs": 0,
        "budget_started_at": None,
        "gave_up": False,
        "last_failure": None,
        "verify_command": None,
    }


def load_state(root: Path, conversation_id: str) -> Dict[str, Any]:
    path = state_path(root, conversation_id)
    if not path.is_file():
        return default_state()
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default_state()
    if not isinstance(loaded, dict):
        return default_state()
    state = default_state()
    state.update(loaded)
    if not isinstance(state.get("files"), dict):
        state["files"] = {}
    return state


def save_state(root: Path, conversation_id: str, state: Mapping[str, Any]) -> None:
    path = state_path(root, conversation_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(state), encoding="utf-8")
    os.replace(temporary, path)


def conversation_id(payload: Mapping[str, Any]) -> str:
    value = payload.get("conversation_id") or payload.get("session_id") or ""
    return str(value)


def watched_rel(root: Path, file_path: str) -> Optional[str]:
    path = Path(file_path)
    if not path.is_absolute():
        path = root / path
    try:
        rel = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return None
    if rel.startswith(".cursor/") or "/vendor/" in "/" + rel or "/build/" in "/" + rel:
        return None
    prefixes = ("crumpling/", "tests/", "notebooks/", "configs/", "legacy/")
    if not rel.startswith(prefixes):
        return None
    if rel.startswith("configs/"):
        return rel
    if rel.endswith((".py", ".ipynb", ".ini", ".cpp", ".h", ".cc", ".hpp")):
        return rel
    return None


def is_moderate(state: Mapping[str, Any]) -> bool:
    files = state.get("files") or {}
    if not isinstance(files, dict) or not files:
        return False
    total = 0
    for value in files.values():
        if isinstance(value, int):
            total += value
    return total >= MODERATE_CHARS or len(files) >= MODERATE_FILES


def handle_edit(payload: Mapping[str, Any]) -> Dict[str, Any]:
    cid = conversation_id(payload)
    if not cid:
        return {}
    root = workspace_root(payload)
    rel = watched_rel(root, str(payload.get("file_path") or ""))
    if rel is None:
        return {}
    chars = 0
    edits = payload.get("edits") or []
    if isinstance(edits, list):
        for edit in edits:
            if not isinstance(edit, dict):
                continue
            old = str(edit.get("old_string") or "")
            new = str(edit.get("new_string") or "")
            chars += max(len(old), len(new))
    if chars == 0:
        chars = 80
    state = load_state(root, cid)
    was_verified = bool(state.get("verified"))
    files = dict(state.get("files") or {})
    files[rel] = int(files.get(rel) or 0) + chars
    state["files"] = files
    state["verified"] = False
    if was_verified:
        state["repairs"] = 0
        state["budget_started_at"] = None
        state["gave_up"] = False
        state["last_failure"] = None
    save_state(root, cid, state)
    return {}


def tool_fields(payload: Mapping[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    tool_input = payload.get("tool_input") or {}
    if isinstance(tool_input, str):
        try:
            tool_input = json.loads(tool_input)
        except json.JSONDecodeError:
            tool_input = {}
    if not isinstance(tool_input, dict):
        tool_input = {}
    command = str(tool_input.get("command") or payload.get("command") or "")
    output = payload.get("tool_output")
    if output is None:
        output = payload.get("output") or ""
    if not isinstance(output, str):
        output = json.dumps(output)
    return command, tool_input, output


def exit_code_from_output(output: str) -> Optional[int]:
    text = output
    for _ in range(2):
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            break
        if isinstance(parsed, dict):
            for key in ("exitCode", "exit_code", "code"):
                value = parsed.get(key)
                if isinstance(value, int):
                    return value
            nested = parsed.get("output") or parsed.get("stdout") or ""
            text = nested if isinstance(nested, str) else json.dumps(parsed)
        elif isinstance(parsed, str):
            text = parsed
        else:
            break
    match = re.search(r"Exit code:\s*(-?\d+)", text)
    if match:
        return int(match.group(1))
    return None


def is_verify_command(command: str) -> bool:
    if "--collect-only" in command or re.search(r"(^|\s)--help(\s|$)", command):
        return False
    patterns = (
        r"\bpytest\b",
        r"\bcrumpling-swapstat\b",
        r"-m\s+crumpling\b",
        r"nbconvert\b[\s\S]{0,240}--execute",
        r"--execute\b[\s\S]{0,240}nbconvert",
        r"\bpapermill\b",
        r"build_core\.py",
    )
    return any(re.search(pattern, command) for pattern in patterns)


def handle_shell(payload: Mapping[str, Any], failed: bool) -> Dict[str, Any]:
    cid = conversation_id(payload)
    if not cid:
        return {}
    command, _, output = tool_fields(payload)
    if not is_verify_command(command):
        return {}
    root = workspace_root(payload)
    state = load_state(root, cid)
    if not is_moderate(state):
        return {}
    code = 1 if failed else exit_code_from_output(output)
    if code == 0:
        state["verified"] = True
        state["verify_command"] = command
        state["last_failure"] = None
    elif code is not None:
        state["verified"] = False
        state["last_failure"] = {
            "command": command,
            "code": code,
            "output": trim(output),
        }
    save_state(root, cid, state)
    return {}


def handle_user_prompt(payload: Mapping[str, Any]) -> Dict[str, Any]:
    prompt = str(payload.get("prompt") or "")
    cid = conversation_id(payload)
    if not cid or FOLLOWUP_MARK in prompt:
        return {"continue": True}
    root = workspace_root(payload)
    state = load_state(root, cid)
    if state.get("budget_started_at") or state.get("gave_up") or state.get("repairs"):
        state["budget_started_at"] = None
        state["gave_up"] = False
        state["repairs"] = 0
        state["last_failure"] = None
        save_state(root, cid, state)
    return {"continue": True}


def handle_stop(payload: Mapping[str, Any]) -> Dict[str, Any]:
    if payload.get("status") == "aborted":
        return {}
    cid = conversation_id(payload)
    if not cid:
        return {}
    root = workspace_root(payload)
    state = load_state(root, cid)
    if not is_moderate(state) or state.get("verified") or state.get("gave_up"):
        return {}
    now = time.time()
    started = state.get("budget_started_at")
    if isinstance(started, (int, float)) and now - float(started) > BUDGET_SECONDS:
        return give_up(root, cid, state, "the 10 minute budget is spent")
    code, command, output = run_checks(root, list(state.get("files") or {}))
    if code == 0:
        state["verified"] = True
        state["verify_command"] = command
        state["last_failure"] = None
        save_state(root, cid, state)
        return {}
    failure = {"command": command, "code": code, "output": trim(output)}
    state["last_failure"] = failure
    if not isinstance(started, (int, float)):
        state["budget_started_at"] = now
    if int(state.get("repairs") or 0) >= REPAIR_LIMIT:
        return give_up(root, cid, state, "3 repair rounds are used")
    state["repairs"] = int(state.get("repairs") or 0) + 1
    save_state(root, cid, state)
    return {"followup_message": repair_message(state, failure)}


def give_up(
    root: Path, cid: str, state: Dict[str, Any], reason: str
) -> Dict[str, Any]:
    state["gave_up"] = True
    save_state(root, cid, state)
    failure = state.get("last_failure") or {}
    command = failure.get("command") or "the relevant check"
    code = failure.get("code")
    output = failure.get("output") or ""
    message = (
        "{mark} Verification stopped because {reason}. "
        "Do not edit more files. Tell the user the command, the error, "
        "and what you already tried.\n\n"
        "Command: {command}\n"
        "Exit code: {code}\n"
        "Changed files:\n{files}\n\n"
        "{output}"
    ).format(
        mark=FOLLOWUP_MARK,
        reason=reason,
        command=command,
        code=code,
        files=format_files(state),
        output=output,
    )
    return {"followup_message": message}


def repair_message(state: Mapping[str, Any], failure: Mapping[str, Any]) -> str:
    round_number = int(state.get("repairs") or 0)
    code = failure.get("code")
    detail = failure.get("output") or "The check did not finish successfully."
    return (
        "{mark} A moderate change is not green yet. "
        "This is repair round {round} of {limit}. "
        "The budget is {minutes} minutes from the first failure, "
        "and each run gets 3 minutes. "
        "Fix the failure and rerun the smallest check that still exercises "
        "the change. Do not ask the user whether to continue. "
        "If this was a timeout, narrow the notebook or simulation instead of "
        "rewriting working code.\n\n"
        "Command: {command}\n"
        "Exit code: {code}\n"
        "Changed files:\n{files}\n\n"
        "{detail}"
    ).format(
        mark=FOLLOWUP_MARK,
        round=round_number,
        limit=REPAIR_LIMIT,
        minutes=BUDGET_SECONDS // 60,
        command=failure.get("command"),
        code=code,
        files=format_files(state),
        detail=detail,
    )


def format_files(state: Mapping[str, Any]) -> str:
    files = state.get("files") or {}
    names = sorted(files)[:15]
    if not names:
        return "(none)"
    return "\n".join("- " + name for name in names)


def needs_pytest(rels: Sequence[str]) -> bool:
    for rel in rels:
        if rel.startswith("configs/"):
            return True
        if rel.endswith(".py") and rel.startswith(("crumpling/", "tests/")):
            return True
    return False


def tests_for(rels: Sequence[str]) -> List[str]:
    selected = []  # type: List[str]

    def add(item: str) -> None:
        if item not in selected:
            selected.append(item)

    for rel in rels:
        if rel.startswith("tests/") and rel.endswith(".py"):
            add(rel)
            continue
        matched = False
        for prefix, test in MODULE_TESTS:
            if rel == prefix or (prefix.endswith("/") and rel.startswith(prefix)):
                add(test)
                matched = True
                break
        if not matched and rel.startswith("crumpling/") and rel.endswith(".py"):
            add("tests")
    if not selected and any(rel.startswith("configs/") for rel in rels):
        add("tests/test_pipeline.py")
    if not selected:
        add("tests")
    return selected


def run_checks(
    root: Path, rels: Sequence[str]
) -> Tuple[Optional[int], str, str]:
    deadline = time.time() + COMMAND_TIMEOUT_SECONDS
    commands = []  # type: List[str]
    outputs = []  # type: List[str]
    if needs_pytest(rels):
        python = python_for("pytest")
        if python is None:
            return (
                1,
                "python -m pytest",
                "Could not find a Python interpreter that imports pytest.",
            )
        cmd = [python, "-m", "pytest", "-q", "--tb=short"] + tests_for(rels)
        code, output = run_cmd(cmd, root, deadline)
        commands.append(shlex.join(cmd))
        outputs.append(output)
        if code != 0:
            return code, "\n".join(commands), "\n".join(outputs)
    notebooks = [rel for rel in rels if rel.endswith(".ipynb")]
    if notebooks:
        python = python_for("nbconvert", prefer="envs/crumpling")
        if python is None:
            return (
                1,
                "python -m jupyter nbconvert --execute",
                "Could not find a Python interpreter that imports nbconvert. "
                "The crumpling conda env is the one that has Jupyter.",
            )
        for notebook in notebooks:
            cmd = [
                python,
                "-m",
                "jupyter",
                "nbconvert",
                "--execute",
                "--to",
                "notebook",
                "--stdout",
                notebook,
            ]
            code, output = run_cmd(cmd, root, deadline)
            commands.append(shlex.join(cmd))
            if code != 0:
                outputs.append(output)
                return code, "\n".join(commands), "\n".join(outputs)
            outputs.append("executed " + notebook)
    if commands:
        return 0, "\n".join(commands), "\n".join(outputs)
    manual = manual_check(rels)
    return None, manual, "Run this yourself within 3 minutes:\n" + manual


def manual_check(rels: Sequence[str]) -> str:
    if any(rel.startswith("legacy/") for rel in rels):
        return (
            "python legacy/cpp_core/build_core.py "
            "or the changed slice of notebooks/compare_cpp_vs_python.ipynb"
        )
    return "python -m pytest -q tests"


def run_cmd(cmd: Sequence[str], root: Path, deadline: float) -> Tuple[int, str]:
    remaining = deadline - time.time()
    if remaining < 5:
        return 124, "The 3 minute check budget was already used."
    try:
        completed = subprocess.run(
            list(cmd),
            cwd=str(root),
            capture_output=True,
            text=True,
            timeout=remaining,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = getattr(exc, "stdout", "") or ""
        stderr = getattr(exc, "stderr", "") or ""
        return 124, trim(str(stdout) + "\n" + str(stderr)) + "\nTimed out after 3 minutes."
    except OSError as exc:
        return 127, str(exc)
    text = completed.stdout or ""
    if completed.stderr:
        text = text + "\n" + completed.stderr
    return completed.returncode, trim(text)


def trim(text: str, limit: int = OUTPUT_LIMIT) -> str:
    if len(text) <= limit:
        return text
    return "...\n" + text[-limit:]


def candidate_pythons() -> List[str]:
    found = []  # type: List[str]
    current = shutil.which("python3")
    if current:
        found.append(current)
    found.append(sys.executable)
    home = Path.home()
    for distro in ("miniconda3", "anaconda3", "mambaforge", "miniforge3"):
        for env_name in ("crumpling", "generalEnv_linux"):
            found.append(
                str(home / distro / "envs" / env_name / "bin" / "python")
            )
    unique = []  # type: List[str]
    for item in found:
        if item not in unique:
            unique.append(item)
    return unique


def python_for(module: str, prefer: str = "") -> Optional[str]:
    cache_path = Path(__file__).resolve().parent / "state" / "interpreters.json"
    cache = {}  # type: Dict[str, str]
    if cache_path.is_file():
        try:
            loaded = json.loads(cache_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            loaded = {}
        if isinstance(loaded, dict):
            cache = {str(key): str(value) for key, value in loaded.items()}
    cache_key = module + "|" + prefer
    cached = cache.get(cache_key)
    if cached and Path(cached).is_file() and probe(cached, module):
        return cached
    candidates = [item for item in candidate_pythons() if Path(item).is_file()]
    if prefer:
        candidates.sort(key=lambda item: 0 if prefer in item else 1)
    for candidate in candidates:
        if probe(candidate, module):
            cache[cache_key] = candidate
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cache), encoding="utf-8")
            return candidate
    return None


def probe(python: str, module: str) -> bool:
    try:
        completed = subprocess.run(
            [python, "-c", "import " + module],
            capture_output=True,
            timeout=20,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return completed.returncode == 0


def handle_guard_write(payload: Mapping[str, Any]) -> Dict[str, Any]:
    root = workspace_root(payload)
    _, tool_input, _ = tool_fields(payload)
    raw_path = tool_input.get("path") or tool_input.get("file_path") or ""
    if not raw_path:
        return {"permission": "allow"}
    target = protected_new_module(root, str(raw_path))
    if target is None:
        return {"permission": "allow"}
    if approved(payload.get("transcript_path"), target):
        return {"permission": "allow"}
    return deny_new_module(target)


def handle_guard_shell(payload: Mapping[str, Any]) -> Dict[str, Any]:
    command = str(payload.get("command") or "")
    if not command or not CREATE_VERB_RE.search(command):
        return {"permission": "allow"}
    root = workspace_root(payload)
    stripped = command.replace('"', "").replace("'", "")
    targets = []  # type: List[str]
    for match in PY_PATH_RE.findall(stripped):
        rel = normalize_repo_path(root, match)
        if rel and protected_new_module(root, rel):
            targets.append(rel)
    if re.search(r"\bmkdir\b", stripped):
        for match in DIR_PATH_RE.findall(stripped):
            if match.endswith(".py"):
                continue
            rel = normalize_repo_path(root, match)
            if rel is None:
                continue
            path = root / rel
            if not path.exists():
                targets.append(rel + "/")
    unique = []  # type: List[str]
    for target in targets:
        if target not in unique:
            unique.append(target)
    if not unique:
        return {"permission": "allow"}
    pending = [
        target
        for target in unique
        if not approved(payload.get("transcript_path"), target.rstrip("/"))
    ]
    if not pending:
        return {"permission": "allow"}
    shown = ", ".join(pending)
    return {
        "permission": "ask",
        "user_message": (
            "This command would create a new crumpling module: {0}. "
            "Approve it only if that new module should exist."
        ).format(shown),
        "agent_message": (
            "Creating {0} needs approval. Prefer an existing module. "
            "If you still need it, wait for the user to approve this command."
        ).format(shown),
    }


def normalize_repo_path(root: Path, raw: str) -> Optional[str]:
    path = Path(raw)
    if path.is_absolute():
        try:
            return path.resolve().relative_to(root.resolve()).as_posix()
        except ValueError:
            text = path.as_posix()
            marker = "/crumpling/"
            index = text.find(marker)
            if index == -1:
                return None
            return text[index + 1 :]
    return path.as_posix()


def protected_new_module(root: Path, raw_path: str) -> Optional[str]:
    path = Path(raw_path)
    if not path.is_absolute():
        path = root / path
    try:
        rel = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return None
    if not rel.startswith("crumpling/") or not rel.endswith(".py"):
        return None
    if path.exists():
        return None
    return rel


def approved(transcript_path: Any, rel: str) -> bool:
    if not transcript_path:
        return False
    path = Path(str(transcript_path))
    messages = load_messages(path)
    proposed_at = None  # type: Optional[int]
    for index, (role, text) in enumerate(messages):
        if role == "assistant" and rel in text and re.search(
            r"creat|new module|approval", text, re.IGNORECASE
        ):
            proposed_at = index
        if role != "user":
            continue
        if FOLLOWUP_MARK in text:
            continue
        body = user_body(text)
        if rel in body and APPROVE_RE.search(body) and not REJECT_RE.search(body):
            return True
        if proposed_at is not None and index > proposed_at:
            stripped = body.strip()
            if (
                len(stripped) <= 80
                and APPROVE_RE.search(stripped)
                and not REJECT_RE.search(stripped)
            ):
                return True
            if REJECT_RE.search(stripped) and len(stripped) <= 80:
                proposed_at = None
    return False


def user_body(text: str) -> str:
    match = re.search(r"<user_query>\s*(.*?)\s*</user_query>", text, re.DOTALL)
    if match:
        return match.group(1)
    return text


def load_messages(path: Path) -> List[Tuple[str, str]]:
    if not path.is_file():
        return []
    try:
        size = path.stat().st_size
        with path.open("rb") as handle:
            if size > 1_500_000:
                handle.seek(size - 1_500_000)
                handle.readline()
            data = handle.read().decode("utf-8", errors="replace")
    except OSError:
        return []
    messages = []  # type: List[Tuple[str, str]]
    for line in data.splitlines():
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(obj, dict):
            continue
        role = str(obj.get("role") or "")
        message = obj.get("message") if isinstance(obj.get("message"), dict) else obj
        content = message.get("content") if isinstance(message, dict) else ""
        text = content_text(content)
        if role and text:
            messages.append((role, text))
    return messages


def content_text(content: Any) -> str:
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        chunks = []  # type: List[str]
        for part in content:
            if isinstance(part, dict) and part.get("type") == "text":
                chunks.append(str(part.get("text") or ""))
        return "\n".join(chunks)
    return ""


def deny_new_module(rel: str) -> Dict[str, Any]:
    homes = ", ".join(EXISTING_MODULES)
    return {
        "permission": "deny",
        "user_message": (
            "Blocked a new module until you approve it: {0}."
        ).format(rel),
        "agent_message": (
            "Do not create {0}. Put the code in an existing module if it fits "
            "({1}). If it truly does not fit, ask the user to approve this exact "
            "path and wait. Retry the write only after they explicitly approve it."
        ).format(rel, homes),
    }


if __name__ == "__main__":
    main()
