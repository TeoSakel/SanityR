# Contributing to SanityR

Thank you for your interest in contributing to SanityR! This document provides
guidelines to help you get started.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Reporting Bugs](#reporting-bugs)

---

## Code of Conduct

This project follows the
[Bioconductor Code of Conduct](https://bioconductor.org/about/code-of-conduct/).
Please be respectful and constructive in all interactions.

---

## How to Contribute

1. **Bug reports** — open an
   [issue](https://github.com/TeoSakel/SanityR/issues) with a minimal
   reproducible example.
2. **Feature requests** — open an issue describing the use-case before
   implementing anything large.
3. **Pull requests** — fork the repository, make your changes on a feature
   branch, and open a PR against `master`.

---

## Development Setup

```r
# Install development dependencies
install.packages(c("devtools", "BiocCheck", "roxygen2"))

# Clone and load the package
git clone https://github.com/TeoSakel/SanityR.git
cd SanityR

# In R
devtools::load_all()
```

After editing C++ files in `src/`, regenerate the Rcpp exports:

```r
Rcpp::compileAttributes()
devtools::load_all()
```

After editing roxygen2 comments, regenerate `man/` and `NAMESPACE`:

```r
devtools::document()
```

---

## Coding Standards

SanityR follows [Bioconductor coding guidelines](https://contributions.bioconductor.org/r-code.html):

- **Indentation**: 4 spaces, no tabs.
- **Line length**: ≤ 80 characters.
- **Assignment**: `<-` (not `=`).
- **Booleans**: `TRUE`/`FALSE` (never `T`/`F`).
- **Sequences**: `seq_len(n)` and `seq_along(x)` (never `1:n`).
- **S4**: use `setGeneric()` / `setMethod()` for new generics.
- **Parallelism**: accept a `BPPARAM` argument; use `bplapply()`.
- **Imports**: add to `DESCRIPTION` only; access via `pkg::fun()` in code.
- **Never** call `library()` or `require()` inside package source files.
- **Never** edit `man/` or `NAMESPACE` by hand — these are auto-generated.

---

## Testing

Tests live in `tests/testthat/` using **testthat edition 3**. Each source file
should have a corresponding `test-<topic>.R` file.

```r
# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-Sanity.R")
```

Guidelines:

- Cover both the happy path and expected errors/warnings.
- For `SingleCellExperiment` outputs, assert assay names and `rowData` columns.
- Use `set.seed()` for any randomised tests to ensure reproducibility.
- Do not add large data files to the test suite; use `simulate_*` helpers or
  the example data in `inst/extdata/`.

---

## Submitting Changes

1. Fork the repo and create a branch: `git checkout -b fix/short-description`.
2. Make your changes following the standards above.
3. Run the full check suite:

   ```r
   devtools::check()        # R CMD check
   BiocCheck::BiocCheck()   # Bioconductor-specific checks
   ```

4. Ensure **no new `ERROR`s, `WARNING`s, or `NOTE`s** are introduced.
5. Open a pull request with a clear description of *what* changed and *why*.

---

## Reporting Bugs

Please include:

- SanityR version (`packageVersion("SanityR")`).
- R version (`R.version.string`).
- Bioconductor version (`BiocManager::version()`).
- A **minimal reproducible example** (preferably using `simulate_*` helpers).
- The full error message or unexpected output.
