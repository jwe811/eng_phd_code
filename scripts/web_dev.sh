#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BACKEND_PID=""

cleanup() {
    if [[ -n "${BACKEND_PID}" ]] && kill -0 "${BACKEND_PID}" 2>/dev/null; then
        kill "${BACKEND_PID}" 2>/dev/null || true
        wait "${BACKEND_PID}" 2>/dev/null || true
    fi
}

trap cleanup EXIT INT TERM

cd "${ROOT}"

if [[ ! -x web/backend/.venv/bin/uvicorn ]]; then
    echo "Backend virtualenv is missing dependencies."
    echo "Run: python3 -m venv web/backend/.venv && web/backend/.venv/bin/pip install -r web/backend/requirements.txt"
    exit 1
fi

if [[ ! -d web/frontend/node_modules ]]; then
    echo "Frontend dependencies are missing."
    echo "Run: cd web/frontend && npm install"
    exit 1
fi

echo "Starting backend at http://127.0.0.1:8000"
PYTHONPATH="${ROOT}:${ROOT}/web/backend" \
    web/backend/.venv/bin/uvicorn app.main:app \
    --app-dir web/backend \
    --host 127.0.0.1 \
    --port 8000 \
    --reload \
    --reload-dir web/backend \
    --reload-dir postprocess &
BACKEND_PID="$!"

echo "Starting frontend at http://127.0.0.1:5173"
cd web/frontend
npm run dev -- --host 127.0.0.1
