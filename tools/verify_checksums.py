#!/usr/bin/env python3
"""Verify files against config/file_checksums_sha256.txt."""

from pathlib import Path
import hashlib
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config/file_checksums_sha256.txt"

failures = []
checked = 0
for line in MANIFEST.read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    digest, rel = line.split("  ", 1)
    path = ROOT / rel
    if not path.exists():
        failures.append(f"missing: {rel}")
        continue
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != digest:
        failures.append(f"checksum mismatch: {rel}")
    checked += 1

if failures:
    print("CHECKSUM VERIFICATION FAILED")
    for item in failures:
        print(" -", item)
    sys.exit(1)

print(f"CHECKSUM VERIFICATION PASSED ({checked} files)")
