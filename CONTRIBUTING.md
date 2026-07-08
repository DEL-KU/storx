# Contributing to STORX

Thank you for your interest in improving STORX. This document explains how to report issues, propose changes, and add new material to the framework.

## Reporting issues

Please open an issue on the [GitHub issue tracker](https://github.com/DEL-KU/storx/issues) and include:

- The MATLAB version and operating system you are using.
- The script or class involved (for example, `05-topopt2d/topopt2d_density.m`).
- Steps to reproduce the problem, and the error message or unexpected output.

## Suggesting educational content

STORX is organized as a sequence of numbered chapters under `00-examples/`. If you would like to propose a new chapter, example, or exercise, please open an issue first describing the topic and how it fits into the existing sequence before submitting a pull request.

## Submitting changes

1. Fork the repository and create a branch for your change.
2. Follow the existing class structure and abstract interfaces (for example, new physics solvers should extend `simulation2d`/`fea2d`, and new optimization methods should follow the patterns used in `05-topopt2d/`).
3. Add or update tests in `tests/` for any new class or method. Run the full suite locally before opening a pull request:
   ```matlab
   run('runMeFirst');
   runtests('tests');
   ```
4. Open a pull request describing the motivation for the change and, if applicable, the chapter or example it supports.

## Code style

- Match the object-oriented conventions already used in the repository (abstract base classes, defined public APIs).
- Keep methods documented with brief comments explaining non-obvious steps, consistent with the existing codebase.
- Prefer readability over micro-optimization, since STORX is primarily an educational framework.

## Questions

For questions that are not appropriate for a public issue, contact amirzend@ku.edu.
