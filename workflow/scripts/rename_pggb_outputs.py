#!/usr/bin/env python3
import sys
from pathlib import Path
import os

def longest_common_prefixA(strings):
    """Return the longest common prefix of a list of strings."""
    if not strings:
        return ""
    # Start with the first string entirely, then chop off until it's a prefix of all.
    strings = [s for s in strings if 'done' not in s]
    strings = [s for s in strings if 'alignments' not in s]

    prefix = strings[0]
    for s in strings[1:]:
        while not s.startswith(prefix) and prefix:
            prefix = prefix[:-1]
    return prefix

def longest_common_prefixB(strings):
    """Return the longest common prefix of a list of strings."""
    if not strings:
        return ""
    # Start with the first string entirely, then chop off until it's a prefix of all.
    strings = [s for s in strings if 'done' not in s]

    prefix = strings[0]
    for s in strings[1:]:
        while not s.startswith(prefix) and prefix:
            prefix = prefix[:-1]
    return prefix

def rename_pggb_outputs(dirpath, new_prefix):
    """
    In dirpath, find all files, determine their common leading text up to the last '.',
    and rename each file by replacing that old prefix with new_prefix + '.'.
    """
    p = Path(dirpath)
    if not p.is_dir():
        print(f"Error: {dirpath!r} is not a directory.", file=sys.stderr)
        sys.exit(1)

    # Gather all files (not directories) in that folder
    files = [f.name for f in p.iterdir() if f.is_file()]
    if not files:
        print(f"No files found in {dirpath!r}. Nothing to do.", file=sys.stderr)
        sys.exit(0)

    # Compute common prefix
    commonA = longest_common_prefixA(files)
    if '.' in commonA:
        # cut back to the last dot (inclusive)
        dot_pos = commonA.rfind('.') + 1
        old_prefixA = commonA[:dot_pos]
    else:
        # no dot at all—just take the whole common
        old_prefixA = commonA

    commonB = longest_common_prefixB(files)
    if '.' in commonB:
        # cut back to the last dot (inclusive)
        dot_pos = commonB.rfind('.') + 1
        old_prefixB = commonB[:dot_pos]
    else:
        # no dot at all—just take the whole common
        old_prefixB = commonB

    if not old_prefixA and not old_prefixB:
        print("Warning: Could not detect a non-empty common prefix. Aborting.", file=sys.stderr)
        sys.exit(1)

    print(f"Detected old prefix: {old_prefixA!r} (A) and {old_prefixB!r} (B)")
    print(f"Renaming to new prefix: {new_prefix!r}")
    print()

    # Rename each file
    for fname in files:
        if not fname.startswith(old_prefixA):
            if not fname.startswith(old_prefixB):
                print(f"Skipping {fname!r} (does not start with {old_prefixA!r} or {old_prefixB!r})")
                continue

            suffix = fname[len(old_prefixB):]          # e.g. "gfa" or "og.lay"
            new_name = f"{new_prefix}.{suffix}"       # e.g. "test_batch_pangenome.gfa"
            old_path = p / fname
            new_path = p / new_name
            print(f"{fname} → {new_name}")
            os.rename(old_path, new_path)

        suffix = fname[len(old_prefixA):]          # e.g. "gfa" or "og.lay"
        new_name = f"{new_prefix}.{suffix}"       # e.g. "test_batch_pangenome.gfa"
        old_path = p / fname
        new_path = p / new_name
        print(f"{fname} → {new_name}")
        try:
            os.rename(old_path, new_path)
        except FileExistsError:
            print(f"Warning: {new_path!r} already exists. Skipping rename.")
        except FileNotFoundError:
            print(f"Warning: {old_path!r} does not exist. Skipping rename.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <dir> <new_prefix>", file=sys.stderr)
        sys.exit(1)

    dirpath, new_prefix = sys.argv[1], sys.argv[2]
    rename_pggb_outputs(dirpath, new_prefix)
