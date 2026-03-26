#!/usr/bin/env python3
"""Switch the active GitHub Copilot Chat model.

Updates the ``github.copilot.chat.languageModel`` value in
``.vscode/settings.json`` using regex so that JSONC comments are
preserved.  VS Code hot-reloads settings.json, so the change takes
effect immediately — no window reload required.

Usage
-----
    python3 tools/set_copilot_model.py <model>

Available <model> aliases
-------------------------
    sonnet   →  claude-sonnet-4-5   (default for this repo)
    opus     →  claude-opus-4-5
    haiku    →  claude-haiku-4-5
    auto     →  (removes the override; Copilot reverts to Auto / GPT-4o)

Examples
--------
    python3 tools/set_copilot_model.py sonnet
    python3 tools/set_copilot_model.py opus
    python3 tools/set_copilot_model.py auto
"""

import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Model aliases → exact Copilot model IDs
# ---------------------------------------------------------------------------
MODELS: dict[str, str] = {
    "sonnet": "claude-sonnet-4-5",
    "opus":   "claude-opus-4-5",
    "haiku":  "claude-haiku-4-5",
}

SETTINGS_FILE = Path(__file__).resolve().parent.parent / ".vscode" / "settings.json"

# Regex that matches the setting with any current value (including JSONC lines)
_PATTERN = re.compile(
    r'("github\.copilot\.chat\.languageModel"\s*:\s*)"[^"]*"'
)


def _read() -> str:
    return SETTINGS_FILE.read_text(encoding="utf-8")


def _write(content: str) -> None:
    SETTINGS_FILE.write_text(content, encoding="utf-8")


def switch_to(model_id: str) -> None:
    """Replace the languageModel value with *model_id*."""
    content = _read()
    new_content, n = _PATTERN.subn(rf'\1"{model_id}"', content)
    if n == 0:
        print(f"Warning: 'github.copilot.chat.languageModel' not found in {SETTINGS_FILE}.")
        print("No changes made.")
        sys.exit(1)
    _write(new_content)
    print(f"✓  Copilot Chat model → {model_id}")
    print("   VS Code picks up the change immediately (no reload needed).")


def remove_override() -> None:
    """Remove the languageModel setting so Copilot falls back to Auto."""
    content = _read()
    # Remove the entire line (handles trailing comma and optional comment on same line)
    new_content = re.sub(
        r'[ \t]*"github\.copilot\.chat\.languageModel"\s*:[^\n]*\n?',
        "",
        content,
    )
    if new_content == content:
        print("No 'github.copilot.chat.languageModel' setting found — already using Auto.")
        return
    _write(new_content)
    print("✓  Removed languageModel override → Copilot will use Auto (GPT-4o).")
    print("   VS Code picks up the change immediately (no reload needed).")


def main() -> None:
    if len(sys.argv) != 2:
        _usage()

    alias = sys.argv[1].lower()

    if alias == "auto":
        remove_override()
    elif alias in MODELS:
        switch_to(MODELS[alias])
    else:
        print(f"Unknown model alias: {alias!r}")
        _usage()


def _usage() -> None:
    aliases = ", ".join(MODELS.keys()) + ", auto"
    print(f"Usage: python3 tools/set_copilot_model.py <model>")
    print(f"Available aliases: {aliases}")
    sys.exit(1)


if __name__ == "__main__":
    main()
