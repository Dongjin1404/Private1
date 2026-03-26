# GitHub Copilot Instructions

This repository is a bioinformatics project (Private1) containing coursework and data-library software for Algo_BI and SS_BI, as well as data-fetching utilities (NCBI, Kaggle).

## How to Check Your Copilot Plan (Pro vs Pro+)

### Quick check — GitHub billing page

1. Go to **[github.com/settings/billing](https://github.com/settings/billing)**.
2. Look for the **GitHub Copilot** line under *Current plan* or *Add-ons*.
   - If it reads **"Copilot Pro"** → you have Copilot Pro.
   - If it reads **"Copilot Pro+"** → you have Copilot Pro+.

### Alternative check — Copilot settings page

1. Go to **[github.com/settings/copilot](https://github.com/settings/copilot)**.
2. Open the **Models** section.  
   - **Copilot Pro**: Only the standard OpenAI models (GPT-4o, etc.) are listed; Claude models are absent or shown as requiring an upgrade.  
   - **Copilot Pro+**: Claude Sonnet 4.5 and Claude Opus 4.5 appear as selectable models (they may be listed but only activated via workspace settings in some UI versions).

### Key differences between the two plans

| Feature | Copilot Pro | Copilot Pro+ |
|---|---|---|
| GPT-4o / GPT-4.1 | ✅ | ✅ |
| Claude Sonnet 4.5 | ❌ | ✅ |
| Claude Opus 4.5 | ❌ | ✅ |
| Gemini models | Limited | ✅ |
| Monthly price (approx.) | $10/mo | $19/mo |

> **Note for this repository:** The workspace is configured to use **Claude Sonnet 4.5** (`"github.copilot.chat.languageModel": "claude-sonnet-4-5"`). This model requires **Copilot Pro+**. If your plan is Copilot Pro, requests will either fall back to the default model or produce an error — upgrade to Copilot Pro+ at [github.com/settings/billing](https://github.com/settings/billing) to use Claude models.

---

## Using Claude Models with Copilot Pro+

Copilot Pro+ includes access to **Claude Sonnet 4.5** and **Claude Opus 4.5** models in addition to the default OpenAI models.  
These are enabled for this workspace. To switch to a Claude model:

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

**Why Claude Sonnet 4.5 / Opus 4.5 do not appear in the model picker dropdown for "Local":**  
The model picker in the web Codespace only lists models that GitHub exposes as selectable in the UI for your current plan.  
Newer model versions (4.6 and above) appear in the list but are greyed out because they require a higher plan tier.  
Claude Sonnet 4.5 and Opus 4.5 are available to Copilot Pro+ but are applied via workspace settings rather than the picker UI.  
The workspace setting `"github.copilot.chat.languageModel": "claude-sonnet-4-5"` in `.devcontainer/devcontainer.json` sets the active model automatically — even when the model picker does not list it explicitly.

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
   - Click the **model picker** at the top of the chat input and select **Claude Opus**, **or**
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
- **Why the model picker shows greyed-out 4.6/5.x models**: Those versions require a plan upgrade. Claude 4.5 models are your plan's included models and are applied via workspace settings, not shown as separate picker entries in every UI version.
- **Why Claude 4.5 does not appear in the picker at all**: GitHub's model picker UI does not always list every available model — it only lists those GitHub has exposed as individually selectable. The workspace setting `"github.copilot.chat.languageModel": "claude-sonnet-4-5"` still routes requests to Claude Sonnet 4.5 even when the name is absent from the dropdown.
- **Web Codespace "Local" defaults to Gemini Pro**: The UI default label may show Gemini Pro or another model name, but the active model is determined by the workspace setting in `.devcontainer/devcontainer.json`. The label updates once you send a message.
- **Local VS Code**: After opening the project, reload the window once (`Ctrl+Shift+P` → "Developer: Reload Window") so VS Code picks up the workspace settings in `.vscode/settings.json`.
- After changing the model setting, reload the Codespace tab in your browser (or run `Ctrl+Shift+P` → "Developer: Reload Window" in the Codespace) for the change to take effect.
- If the model picker does not show Claude models, sign out and back into GitHub in VS Code / the web Codespace:
  `Ctrl+Shift+P` → "GitHub: Sign Out", then "GitHub: Sign In".
