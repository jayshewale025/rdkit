#!/usr/bin/env bash
# CPython 3.13 ARM64 wheel (Raspberry Pi). Kept for the existing workflow; see build_linux.sh.
set -euo pipefail

test "$(uname -m)" = aarch64
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY_TAG=cp313 BUILD_DIR="$root/build-mechanism-arm" JOBS="${JOBS:-2}" \
  exec bash "$root/Packaging/mechanism/build_linux.sh"
