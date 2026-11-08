#!/usr/bin/env bash
# Fetch annexed MRI files over public HTTPS (OpenNeuro S3 object URLs recorded
# in the annex), verify SHA-256 against the git-annex key, then reinject.
# Used because git-annex's S3 client is rejected by this sandbox's proxy
# credential injection. Usage: scripts/fetch_mri.sh <path> [<path> ...]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
tmp="${TMPDIR:-/tmp}/annex_fetch.$$"; mkdir -p "$tmp"; trap 'rm -rf "$tmp"' EXIT
for f in "$@"; do
  [ -e "$f" ] && { echo "present  $f"; continue; }
  url=$(git annex whereis --json "$f" | python3 -c 'import json,sys;d=json.load(sys.stdin);print([u for w in d["whereis"] for u in w.get("urls",[]) if u.startswith("http")][0])')
  key=$(basename "$(readlink "$f")")
  want=$(echo "$key" | sed -E 's/^SHA256E-s[0-9]+--([0-9a-f]{64}).*/\1/')
  out="$tmp/$key"
  for i in 1 2 3 4; do curl -sSfL --retry 3 -o "$out" "$url" && break; sleep $((2**i)); done
  got=$(sha256sum "$out" | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "CHECKSUM MISMATCH $f" >&2; exit 1; }
  git annex reinject "$out" "$f" >/dev/null
  echo "ok sha256=$got  $f"
done
