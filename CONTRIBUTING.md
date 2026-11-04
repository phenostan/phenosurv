# Contributing to phenosurv

Thanks for helping build phenosurv! This guide covers the development environment, how we write code here (tests first), and how to review it. For what the package does and how to use it, start with the [README](README.md).

- [Local setup](#local-setup)
- [Test-driven development](#test-driven-development)
- [Reviewing code](#reviewing-code)
- [Submitting changes](#submitting-changes)

## Local setup

We use [pixi](https://pixi.sh) to manage the development environment. One command installs Python, [CmdStan](https://mc-stan.org/docs/cmdstan-guide/) (including the C++ compiler it needs) and every other dependency, at the exact versions recorded in `pixi.lock`. You don't need conda, a virtual environment or `cmdstanpy.install_cmdstan()`.

1. **Install pixi** (once per computer):

   - macOS / Linux:
     ```sh
     curl -fsSL https://pixi.sh/install.sh | sh
     ```
   - Windows (PowerShell):
     ```powershell
     powershell -ExecutionPolicy ByPass -c "irm -useb https://pixi.sh/install.ps1 | iex"
     ```

   Restart your terminal afterwards so the `pixi` command is available.

2. **Clone the repository and install the environment:**
   ```sh
   git clone https://github.com/phenostan/phenosurv.git
   cd phenosurv
   pixi install
   ```
   The first install downloads CmdStan and a compiler, so it can take a few minutes. The environment lives in `.pixi/` inside the repo and doesn't affect anything else on your computer.

3. **Check that everything works:**
   ```sh
   pixi run test      # unit tests
   pixi run example   # end-to-end fit, diagnose and predict (examples/worked_example.py)
   ```

### Day-to-day commands

| What | Command |
| --- | --- |
| Run all tests | `pixi run test` |
| Run one test file | `pixi run python -m unittest tests.test_thresh_surv_gdd -v` |
| Run one test | `pixi run python -m unittest tests.test_thresh_surv_gdd.TestThreshSurvGDDPrintCode.test_print_code_writes_to_stdout` |
| Run the worked example and regenerate `docs/figures/` | `pixi run example` |
| Regenerate the simulated test data | `pixi run simulate-testdata` (see `--help` for options) |
| Update the packaged Stan models from phenostan | `pixi run sync-models` (add `--ref v0.1.0` for a specific tag, branch or commit) |
| Open a shell inside the environment | `pixi shell` (then use `python` directly; `exit` to leave) |

The Stan models in `src/phenosurv/stan/` are copied from [phenostan](https://github.com/phenostan/phenostan); don't edit them here. Change the model in phenostan, then run `pixi run sync-models` and commit the updated files together with `PHENOSTAN_VERSION`.

phenosurv itself is installed in editable mode, so changes under `src/phenosurv/` take effect right away without reinstalling.

### Adding a dependency

- `pixi add <package>` adds a conda-forge package. Prefer this, especially for anything that needs compiled code.
- `pixi add --pypi <package>` adds a package that is only on PyPI.
- If phenosurv needs the package **at runtime** (not just for development), also add it to `dependencies` under `[project]` in `pyproject.toml`. That list is what people get when they `pip install phenosurv`.
- Commit the updated `pyproject.toml` **and** `pixi.lock` together.

## Test-driven development

phenosurv is built test-first. Every change to behavior starts with a test that fails, and code is only written to make a failing test pass. We use the standard-library [`unittest`](https://docs.python.org/3/library/unittest.html) framework.

### The cycle: red → green → refactor

1. **Red.** Write a small test for the next behavior you want. Run it and **watch it fail**, and check that it fails for the reason you expect (for example `ImportError: cannot import name 'ThreshSurvGDD'`, not a typo). A test you've never seen fail might not be testing anything.
2. **Green.** Write the simplest code that makes the test pass. Resist adding features that no test asks for yet.
3. **Refactor.** With all tests passing, clean up names, duplication and structure. Rerun the tests after each change.

Then pick the next behavior and repeat.

### What makes a good test here

- **Test behavior, not implementation.** Check what a user of phenosurv sees: return values, printed output, raised errors. Don't assert on how the code calls cmdstanpy internally, so the implementation can change without rewriting tests.
- **One behavior per test**, with a name that says what it checks: `test_can_instantiate_without_data` beats `test_model_2`. When a test fails, its name should tell you what broke.
- **Keep tests fast.** Compiling a Stan model or running MCMC takes seconds to minutes. Most tests should check behavior that doesn't need either, so the whole suite stays quick enough to run after every change.
- **Bug fixes start with a test.** Before fixing a bug, write a test that reproduces it and fails. The fix is done when that test passes, and the test keeps the bug from coming back.

### Where tests go

- Library code lives in `src/phenosurv/`; tests live in `tests/`.
- Test files are named `test_<thing>.py` (for example `tests/test_thresh_surv_gdd.py`) so `unittest` discovers them automatically.
- Group related tests in `unittest.TestCase` classes, and use `setUp` for shared setup.

## Reviewing code

You don't need to know Stan or Bayesian statistics to review most of phenosurv: the Python layer is a thin wrapper around [cmdstanpy](https://mc-stan.org/cmdstanpy/), and the README has a short glossary of the terms used. The Stan models themselves are reviewed in [phenostan](https://github.com/phenostan/phenostan), not here.

**Get the code running first.** Check out the branch, then run `pixi install`, `pixi run test` and `pixi run example`. If any of these fail on your machine, that is itself a useful finding.

**Where the code is.** Almost everything is in `src/phenosurv/thresh_surv_gdd.py` (the `ThreshSurvGDD` class) and its tests in `tests/test_thresh_surv_gdd.py`. `scripts/` holds developer tools that are not part of the installed package.

**What to look for:**

- **Correctness.** Do the methods do what their docstrings say? Are the data passed to Stan in the right shape and with 1-based indices?
- **Tests.** Is each behavior covered by a test that would fail if the behavior broke? Are tests checking behavior rather than implementation (see above)? Is anything important untested?
- **Usability for ecologists.** Are the names, error messages and return types clear to someone who knows their data but not Stan? Would a mistake in the input data be caught early with a helpful message?
- **Packaging.** Does the package work when installed rather than run from the repo (e.g. Stan files found through `importlib.resources`, compiled executables written somewhere writable)?
- **Readability.** Naming, docstrings, duplication and anything that made the code hard to follow.

When you leave a comment, say whether it must be fixed before merging or is a suggestion, and point to the file and line.

## Submitting changes

1. Create a branch for your work: `git switch -c short-description`.
2. Work through red → green → refactor, committing as you go.
3. Make sure `pixi run test` passes.
4. Open a pull request against `main`. In the description, link the issue it addresses (for example `Closes #2`) and summarize the behavior the new tests cover.
