#!/bin/bash
# Claude Code on the web のセッション開始時に依存関係を用意する。
# ローカル（VS Code など）では何もしない。ローカルは .venv を各自で作る想定。
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}"

PYTHON="$(command -v python3 || command -v python)"

# リモートコンテナは使い捨てなので venv は作らず、そのまま入れる。
"$PYTHON" -m pip install --quiet --disable-pip-version-check \
  --root-user-action=ignore -r requirements-dev.txt

# cce_bg.py をリポジトリ直下から import できるようにする。
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo "export PYTHONPATH=\"${CLAUDE_PROJECT_DIR:-$PWD}\"" >> "$CLAUDE_ENV_FILE"
fi

echo "session-start: pytest / ruff / numpy の準備が完了しました"
