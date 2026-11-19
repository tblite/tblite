#!/usr/bin/env python3
"""Pass when a command rejects its input with a normal error exit."""

import subprocess
import sys


command = sys.argv[1:]
if not command:
    raise SystemExit("expected a command to run")

result = subprocess.run(command, capture_output=True, text=True)
output = result.stdout + result.stderr

if result.returncode != 1:
    sys.stderr.write(output)
    raise SystemExit(f"expected exit status 1, got {result.returncode}")
if not output.strip():
    raise SystemExit("command failed without reporting an error")

print("Command rejected invalid input with exit status 1")
