---
id: build-and-ci
name: Build, CMake and CI
paths:
  - "**/CMakeLists.txt"
  - "res/cmake/**"
  - "res/ci/**"
  - ".github/**"
---

- A test registered at a new level or with a new label changes where it runs (PR,
  nightly, weekly); say so.
- A new test binary or template instantiation adds compile time to every PR build;
  estimate it.
- Workflows run on pull requests from forks: never expose secrets or write tokens
  to code from the PR.
