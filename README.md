# Option Pricing Lab — C++17

A compact numerical finance project developed by refactoring my Baruch C++ option-pricing coursework. It compares analytical pricing with Monte Carlo simulation and a finite-difference PDE solver, and validates sensitivities and perpetual American exercise boundaries.

The entire implementation is in `option_pricing.cpp`. No Boost, Excel integration, external data, or platform-specific paths are required.

## What the project demonstrates

- Generalized Black–Scholes European call/put prices, analytical delta and gamma.
- Exact geometric Brownian motion terminal simulation with antithetic variates.
- Online variance estimation, standard errors and approximate 95% confidence intervals.
- An explicit finite-difference solver with upwind drift and automatically selected time steps.
- Perpetual American prices with both continuation and immediate-exercise regions.
- Deterministic analytical checks, statistical benchmarks and grid refinement checks.
- CSV output for reproducible comparisons across spot prices and option types.

## Build and run

Requires a C++17 compiler. From the repository root:

```sh
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic option_pricing.cpp -o option_pricing
./option_pricing --test
./option_pricing --demo
./option_pricing --demo > results.csv
```

Or use CMake 3.16 or newer:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The GitHub Actions workflow is configured to build and validate on Linux, macOS and Windows. Those hosted runs will execute after the repository is uploaded; local validation used Apple Clang on macOS.

## Model and methodology

Under the risk-neutral measure, the underlying follows `dS = b S dt + sigma S dW`. The risk-free discount rate is `r`, and `b` is the cost of carry. For a stock with continuous dividend yield `q`, `b = r - q`. Rates and volatility are annualized decimal inputs, and maturity is measured in years.

European pricing uses the generalized Black–Scholes formula. The parity relation is `C - P = S exp((b-r)T) - K exp(-rT)`; this distinction matters when carry differs from the interest rate.

Monte Carlo uses the exact terminal distribution, avoiding Euler time-discretization bias for European payoffs. The demo uses 100,000 independent antithetic pairs (200,000 terminal draws) and seed 42. Each pair average is one observation in Welford's variance estimator. `mc_se` therefore measures uncertainty in the average of independent pairs, not in 200,000 independent paths. The displayed interval uses `estimate ± 1.96 × SE` and is an asymptotic sampling interval.

The PDE is `V_tau = 0.5 sigma² S² V_SS + b S V_S - r V`. It uses a uniform spot grid, an explicit time scheme, and drift stencils with nonnegative off-diagonal weights. Time steps enforce a nonnegative central weight. This improves robustness but introduces first-order spatial drift error. The spot domain ends at `4 × max(S, K)`; call and put boundaries use large-spot asymptotics. Prices at the requested spot use linear interpolation.

The perpetual American engine assumes `r > 0`, `sigma > 0` and `b < r`. Maturity is ignored by this engine. It solves the characteristic quadratic and applies intrinsic value beyond the optimal exercise boundary.

## Validation and results

`--test` checks:

1. The standard call benchmark: `S=K=100`, `r=b=0.05`, `sigma=0.2`, `T=1` gives `10.4505835721856`.
2. Generalized put–call parity and finite-difference Greeks, including zero and negative carry.
3. Monte Carlo agreement within six estimated standard errors plus a small numerical tolerance.
4. PDE agreement within 0.10 price units for the selected benchmarks, and improvement under grid refinement.
5. Expiry, deterministic zero-volatility pricing, invalid input rejection, and undefined Greeks at the payoff kink.
6. Perpetual American reference pricing and call/put immediate exercise regions.

`sample_results.csv` contains the local demo output. `VALIDATION.md` records the actual local checks. The Monte Carlo benchmark is a smoke check, not proof of coverage or convergence for every parameter set. Standard-library normal distributions may produce different draws across platforms despite the same seed.

## Design choices and limitations

The single-file layout makes the work easy to inspect and compile as an application code sample. Engines use value types, scoped random generators and vectors instead of global pointers or manual allocation. Demo inputs can be edited in `quant::demo()`; the command-line interface selects demo, validation or help modes.

This is an educational numerical implementation. It does not include calibration, market data, discrete dividends, transaction costs, finite-maturity American pricing, or exotic/path-dependent payoffs. Large maturities or volatilities can require a wider PDE domain. The time-step rule controls stencil monotonicity; it does not eliminate truncation error. The PDE rejects workloads above its configured budget. Extreme finite parameters can still overflow exponentials. Analytical delta/gamma are reported as NaN at the deterministic payoff kink.

## Coursework provenance and contribution statement

The source material was a Baruch C++ coursework folder organized into Groups A–F, with analytical European and perpetual American pricing, Monte Carlo experiments, and finite-difference exercises. It also included course/example utilities; one supplied Monte Carlo file identifies Datasim Education BC (2008–2011).

This portfolio version reorganizes the numerical concepts into a standalone implementation and adds validation, uncertainty reporting, explicit exercise-region handling, portable builds and automated checks. The coursework folder, supplied utility code, assignment documents, spreadsheets and binaries are not bundled in this repository. The original coursework and third-party examples should retain their own attribution.

This version was refactored with AI assistance. It should be presented as a coursework-derived portfolio project, rather than an entirely unaided implementation. Before using it in an application, review the implementation and be prepared to explain the pricing formulas, antithetic standard error, PDE stencil and American exercise boundaries. If the original work was collaborative, describe your individual contribution accurately.
