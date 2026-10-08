# Local validation

Run date: 2026-10-08 (Asia/Shanghai).
Environment: macOS, Apple Clang 17.0.0, C++17.

- Direct optimized build with `-Wall -Wextra -Wpedantic`: passed with no warnings.
- Built-in analytical, statistical, PDE refinement and boundary checks: passed.
- CMake Release build and CTest: 1/1 tests passed.
- AddressSanitizer + UndefinedBehaviorSanitizer build, running `--test`: passed with no reported diagnostics.
- Demo generated six call/put comparison rows in `sample_results.csv`.

For the at-the-money dividend-paying call in the demo, the analytical value is 9.227006, the Monte Carlo estimate is 9.240241 with SE 0.023051, and the 300-cell PDE value is 9.260974. Numerical approximations are intentionally reported alongside their analytical benchmark.

GitHub-hosted and Windows/Linux builds have not been executed locally; the workflow will validate those environments after upload.
