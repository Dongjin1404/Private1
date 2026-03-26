# Private1

Bioinformatics coursework and data-library software for Algo_BI and SS_BI, plus data-fetching utilities (NCBI, Kaggle).

## Copilot Pro+ — Using Claude Sonnet 4.5 / Opus 4.5

This repository is pre-configured so that **Claude Sonnet 4.5** is the active Copilot Chat model as soon as you open the project — no manual picker selection needed.

### Quick-start

| Where you work | What to do |
|---|---|
| **Browser Codespace** (github.com → Code → Codespaces) | Open the Codespace; select **Local** in the Copilot Chat context picker. Claude Sonnet 4.5 is active automatically. |
| **Local VS Code** | Open the folder; sign in to Copilot Pro+. Claude Sonnet 4.5 is active automatically. |
| **Copilot CLI** (`gh copilot`) | Set your preferred model at [github.com/settings/copilot](https://github.com/settings/copilot). |

### What's included with Copilot Pro+

| Model | Included in Copilot Pro+? |
|---|---|
| Claude Haiku 4.5 | ✅ Yes |
| **Claude Sonnet 4.5** | ✅ Yes |
| **Claude Opus 4.5** | ✅ Yes |
| Claude Sonnet 4.6 | ❌ No — requires Copilot Enterprise/Business |
| Claude Opus 4.6 | ❌ No — requires Copilot Enterprise/Business |

> **Note:** Claude 4.6+ models appear in the model picker as greyed-out upgrade prompts. They are **not** part of Copilot Pro+ — the greyed-out entries are upgrade advertisements, not plan features.

### Why Sonnet 4.5 / Opus 4.5 don't show in the model picker

The VS Code model picker lists models that GitHub has exposed as individually selectable in that UI version.  
In many UI versions the "Claude" section shows only **Haiku 4.5**; Sonnet 4.5 and Opus 4.5 are omitted from the list even though they are fully available to Copilot Pro+ subscribers.

**The workspace `languageModel` setting bypasses the picker entirely.**  
Both `.devcontainer/devcontainer.json` (Codespaces) and `.vscode/settings.json` (local VS Code) already contain:

```json
"github.copilot.chat.languageModel": "claude-sonnet-4-5"
```

This routes every Copilot Chat request to Claude Sonnet 4.5 regardless of what the picker displays.

### Switching models for specific tasks

Because Sonnet 4.5 / Opus 4.5 do not appear in the picker, this repository includes **VS Code Tasks** and a helper script so you can switch models with a single command — no manual file editing required.

#### Option A — VS Code Task (Command Palette, works in browser Codespace and local VS Code)

1. Open the **Command Palette**: `Ctrl+Shift+P` (Windows/Linux) or `Cmd+Shift+P` (macOS)
2. Type **"Tasks: Run Task"** and press Enter
3. Choose one of:
   - **Copilot: Use Claude Sonnet 4.5** — fast and capable (repo default)
   - **Copilot: Use Claude Opus 4.5** — most capable, best for complex reasoning
   - **Copilot: Use Claude Haiku 4.5** — fastest, good for quick completions
   - **Copilot: Use Auto** — reverts to Copilot Auto mode

The task updates `.vscode/settings.json` and VS Code picks up the change **immediately** — no reload required.

#### Option B — Terminal (integrated terminal or Codespace terminal)

```bash
# Switch to Opus for a complex task
python3 tools/set_copilot_model.py opus

# Switch back to Sonnet when done
python3 tools/set_copilot_model.py sonnet

# Available aliases: sonnet, opus, haiku, auto
```

#### Option C — Edit `.vscode/settings.json` directly

Change the one-line value and save:

```json
"github.copilot.chat.languageModel": "claude-opus-4-5"
```

### Detailed instructions and troubleshooting

See [`.github/copilot-instructions.md`](.github/copilot-instructions.md) for:
- Full step-by-step setup for each environment
- How to verify which model is actually active
- Troubleshooting the "Haiku-only" picker, greyed-out 4.6 entries, and plan confusion
