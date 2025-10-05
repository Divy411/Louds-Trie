
## Changes
- The original `main.cpp` used random key generation for testing.
-I replaced it with a CLI program (`merge_cli.cpp`) that accepts two sorted key files and merges them.

## Requirements
- `g++` with C++17 support
-  ensure `-mbmi2` is supported by your CPU

## Build
Clean and rebuild the project:

```bash
make clean
make


./merge_cli keys.txt keys2.txt