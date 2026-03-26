# GitHub Copilot Instructions

This repository is a bioinformatics project (Private1) containing coursework and data-library software for Algo_BI and SS_BI, as well as data-fetching utilities (NCBI, Kaggle).

## Using Claude Models with Copilot Pro+

### What Copilot Pro+ actually includes

| Model | Copilot Pro+ | Higher plan required |
|---|---|---|
| Claude Haiku 4.5 | ✅ Included | — |
| **Claude Sonnet 4.5** | ✅ Included | — |
| **Claude Opus 4.5** | ✅ Included | — |
| Claude Sonnet 4.6 | ❌ Not included | Copilot Enterprise / Business |
| Claude Opus 4.6 | ❌ Not included | Copilot Enterprise / Business |

> **Important:** Copilot Pro+ gives you access to Claude **4.5** models (Haiku, Sonnet, and Opus).
> Claude **4.6 and above** are shown as greyed-out upgrade prompts in the picker and require a higher plan tier — they are **not** included with Copilot Pro+.

### Why Sonnet 4.5 / Opus 4.5 don't appear in the model picker

The VS Code model picker UI has two separate sections:

- A **"Claude" section** that lists models GitHub has explicitly exposed as individually selectable in that UI version — currently this typically shows only **Haiku 4.5** as the directly selectable Claude model.
- A **greyed-out "upgrade" section** that shows 4.6+ models requiring a plan upgrade.

**Neither of these sections shows Sonnet 4.5 or Opus 4.5 as named entries** in many UI versions, even though both are fully available to Copilot Pro+ subscribers. This is a GitHub UI limitation, not a plan restriction.

The correct way to use Sonnet 4.5 or Opus 4.5 is via the **workspace `languageModel` setting**, which bypasses the picker entirely and routes all Copilot Chat requests to the specified model. This repository is already configured to do this automatically.

### How this repository is configured

- `.devcontainer/devcontainer.json` — sets `"github.copilot.chat.languageModel": "claude-sonnet-4-5"` for web Codespaces
- `.vscode/settings.json` — sets the same for local VS Code

When you open this project (locally or in a Codespace), Claude Sonnet 4.5 is active automatically. You do not need to select anything in the model picker.

---

### In a Codespace opened in the browser (github.com → Code → Codespaces)

> **This is the recommended path if you do not have VS Code installed locally.**

The `.devcontainer/devcontainer.json` in this repository already sets **Claude Sonnet 4.5** as the default model. When you open the Codespace the model is applied automatically — you do not need to change anything.

#### Understanding the chat context picker ("Local / Copilot CLI / Cloud")

When you open the Copilot Chat panel in the web Codespace, you may see a dropdown that lets you pick between **Local**, **Copilot CLI**, and **Cloud** chat contexts.

| Context | What it uses |
|---|---|
| **Local** | The VS Code extension running inside the Codespace container. Uses the model set in `settings.json` / `devcontainer.json` — **Claude Sonnet 4.5** for this repo. |
| **Copilot CLI** | The `gh copilot` CLI tool. Uses the model set in your GitHub account settings at [github.com/settings/copilot](https://github.com/settings/copilot). |
| **Cloud** | GitHub Copilot's cloud chat. Uses your account's default model. |

> **Select "Local"** to use the workspace-configured Claude Sonnet 4.5 model. The model picker inside "Local" mode may still show Haiku 4.5 or another model as the "selected" item, but the workspace setting takes precedence — your requests are actually handled by Claude Sonnet 4.5.

#### How to verify Claude Sonnet 4.5 is active in "Local" chat

1. Open the Copilot Chat panel in the browser Codespace.
2. Select **Local** from the context picker.
3. Type `/help` or ask any question. The response header or the model label near the input box should reflect the configured model.
4. Alternatively, open the **Command Palette** (`F1` or `Ctrl+Shift+P`), type **"Open User Settings (JSON)"** and check that `"github.copilot.chat.languageModel"` is `"claude-sonnet-4-5"`.

#### Switching to Claude Opus 4.5 in the web Codespace

The model picker may not list Claude Opus 4.5 directly. To switch:

1. In the Codespace file explorer, open `.vscode/settings.json`.
2. Change:
   ```json
   "github.copilot.chat.languageModel": "claude-sonnet-4-5"
   ```
   to:
   ```json
   "github.copilot.chat.languageModel": "claude-opus-4-5"
   ```
3. Save the file. The change takes effect immediately in the current Codespace session.

### In local VS Code (recommended if you have VS Code installed)
The `.vscode/settings.json` in this repository sets Claude Sonnet 4.5 as the default model automatically when you open the project locally. No extra setup is needed beyond signing in to GitHub Copilot.

1. Open the project folder in VS Code on your local machine.
2. Sign in to GitHub with your Copilot Pro+ account:
   `Ctrl+Shift+P` → **GitHub: Sign In**.
3. Open the Copilot Chat panel (`Ctrl+Alt+I` / `Cmd+Alt+I`).
4. The default model is already set to **Claude Sonnet 4.5** (`claude-sonnet-4-5`).
5. To switch to **Claude Opus 4.5** (`claude-opus-4-5`):
   - Edit `.vscode/settings.json` and change the value of `"github.copilot.chat.languageModel"` to `"claude-opus-4-5"`.

### In the integrated terminal (CLI agent / `gh copilot`)
1. Make sure the GitHub CLI is installed and authenticated:
   ```bash
   gh auth status
   ```
2. Use `gh copilot suggest` or `gh copilot explain` — model selection follows your
   account's default model set in [github.com/settings/copilot](https://github.com/settings/copilot).
3. To use **Claude Sonnet 4.5** or **Claude Opus 4.5** in the CLI, visit:
   **github.com → Settings → Copilot → Models** and set Claude Sonnet 4.5 or Opus 4.5 as your preferred model.

### Troubleshooting
- Make sure your Copilot Pro+ subscription is active at [github.com/settings/billing](https://github.com/settings/billing).
- Claude Sonnet 4.5 and Opus 4.5 require **Copilot Pro+** (not Copilot Pro or Copilot Free).
- **Why the model picker's "Claude" section shows only Haiku 4.5**: The UI exposes only Haiku 4.5 as a directly selectable Claude model in many versions. Sonnet 4.5 and Opus 4.5 are still fully available — use the `languageModel` workspace setting (already configured in this repo) to activate them.
- **Why the model picker shows greyed-out 4.6/5.x models**: Those versions require a plan upgrade beyond Copilot Pro+. They are not included in your plan. Claude 4.5 models are your plan's top-tier included Claude models.
- **Why Copilot Pro+ does NOT include Claude 4.6**: Claude 4.6 and above are available only on Copilot Enterprise or Business plans. Copilot Pro+ stops at the 4.5 generation for Claude models.
- **Why Claude 4.5 does not appear in the picker at all**: GitHub's model picker UI does not always list every available model — it only lists those GitHub has exposed as individually selectable. The workspace setting `"github.copilot.chat.languageModel": "claude-sonnet-4-5"` still routes requests to Claude Sonnet 4.5 even when the name is absent from the dropdown.
- **Web Codespace "Local" defaults to Gemini Pro or another model in the picker UI**: The UI default label may show a different model name, but the active model is determined by the workspace setting in `.devcontainer/devcontainer.json`. The actual model used for your requests is Claude Sonnet 4.5.
- **Local VS Code**: After opening the project, reload the window once (`Ctrl+Shift+P` → "Developer: Reload Window") so VS Code picks up the workspace settings in `.vscode/settings.json`.
- After changing the model setting, reload the Codespace tab in your browser (or run `Ctrl+Shift+P` → "Developer: Reload Window" in the Codespace) for the change to take effect.
- If the model picker does not show Claude models, sign out and back into GitHub in VS Code / the web Codespace:
  `Ctrl+Shift+P` → "GitHub: Sign Out", then "GitHub: Sign In".
