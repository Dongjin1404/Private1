# GitHub Copilot Instructions

This repository is a bioinformatics project (Private1) containing coursework and data-library software for Algo_BI and SS_BI, as well as data-fetching utilities (NCBI, Kaggle).

## Using Claude Models with Copilot Pro+

Copilot Pro+ includes access to **Claude Sonnet** and **Claude Opus** models in addition to the default OpenAI models.  
These are enabled for this workspace. To switch to a Claude model:

### In VS Code / Codespaces (Chat panel)
1. Open the Copilot Chat panel (`Ctrl+Alt+I` / `Cmd+Alt+I`).
2. Click the **model picker** (the model name shown at the top of the chat input, e.g. "GPT-4o").
3. Select **Claude Sonnet** or **Claude Opus** from the dropdown.

### In the integrated terminal (CLI agent / `gh copilot`)
1. Make sure the GitHub CLI is installed and authenticated:
   ```bash
   gh auth status
   ```
2. Use `gh copilot suggest` or `gh copilot explain` — model selection follows your
   account's default model set in [github.com/settings/copilot](https://github.com/settings/copilot).
3. To change the default model for CLI interactions, visit:
   **github.com → Settings → Copilot → Models** and choose Claude Sonnet or Opus as your preferred model.

### Troubleshooting
- Make sure your Copilot Pro+ subscription is active at [github.com/settings/billing](https://github.com/settings/billing).
- Claude models require **Copilot Pro+** (not Copilot Pro or Copilot Free).
- After changing the model in GitHub account settings, reload the Codespace (`Ctrl+Shift+P` → "Developer: Reload Window").
- If the model picker does not show Claude models, sign out and back into GitHub in VS Code:
  `Ctrl+Shift+P` → "GitHub: Sign Out", then "GitHub: Sign In".
