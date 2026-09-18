#!/usr/bin/env bash
# Remove accumulated cruft from biokgsuite-bib.
#
# The agent sandbox cannot delete on this mount, so everything it retired was
# moved to _to_delete/ instead. This script does the actual removal, plus the
# caches and OS files that regenerate on their own.
#
#   bash tools/cleanup.sh          # show what would go, delete nothing
#   bash tools/cleanup.sh --yes    # actually delete
set -euo pipefail
cd "$(dirname "$0")/.."
DRY=1; [ "${1:-}" = "--yes" ] && DRY=0

say() { printf '%s\n' "$*"; }
go()  { if [ "$DRY" = 1 ]; then say "  would remove: $1"; else rm -rf -- "$1"; say "  removed: $1"; fi; }

say "== retired work (_to_delete), $(du -sh _to_delete 2>/dev/null | cut -f1)"
[ -d _to_delete ] && go _to_delete

say "== python caches (regenerate on next run)"
while IFS= read -r d; do go "$d"; done < <(find . -type d -name __pycache__ 2>/dev/null)
[ -d .pytest_cache ] && go .pytest_cache

say "== macOS cruft"
while IFS= read -r f; do go "$f"; done < <(find . -name .DS_Store 2>/dev/null)

say
if [ "$DRY" = 1 ]; then
  say "Dry run. Re-run with --yes to delete."
else
  say "Done. Repo is now $(du -sh . | cut -f1)."
fi
