#!/usr/bin/env python3
"""Generate SHA-256 checksums for the release package."""

from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "config/file_checksums_sha256.txt"

excluded = {
    OUT.resolve(),
}
rows = []
for path in sorted(p for p in ROOT.rglob("*") if p.is_file()):
    if path.resolve() in excluded:
        continue
    rel = path.relative_to(ROOT).as_posix()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows.append(f"{digest}  {rel}")

OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
print(f"Wrote {len(rows)} checksums to {OUT}")
