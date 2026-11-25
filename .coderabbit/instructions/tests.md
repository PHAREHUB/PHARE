---
id: tests
name: Tests that prove the claim
paths:
  - "tests/**"
  - "pyphare/pyphare_tests/**"
---

For each test the PR adds or changes:
- for a bug fix or a new feature: would the test fail if the bug came back or
  the feature were removed? If not, say which assertion is too weak. A test that
  protects existing behaviour during a refactor is expected to pass on master;
- does it check the claim the PR makes, at the precision the PR claims? A
  tolerance much looser than the claimed precision is a finding;
- comparisons must fail on NaN. These pass silently on NaN:
  `np.testing.assert_array_equal` and `assert_allclose` when both sides hold NaN
  at the same places (both treat matching NaNs as equal), `np.allclose(..., equal_nan=True)`,
  and inverted checks such as `assert not x > tol` or `if err > tol: fail()`.
  These fail on NaN, as wanted: `assert x < tol`, `np.allclose` by default,
  googletest `EXPECT_NEAR`;
- floating-point results compared bit for bit across runs are only reproducible
  with 1 or 2 MPI ranks (issue #1526); with 3 or more, compare them with a
  tolerance. Integer results (counts, indices, boxes) can still be compared exactly;
- runtime: a test longer than about 10 minutes on the slowest PR build belongs in
  the nightly tier. Python tests: `phare_test_labels(<id> nightly)` in CMake
  (it only matches `py3_<id>` names). C++ tests have no helper: set the `LABELS`
  property on the test named after the binary, inside `if(TEST <binary>)`
  because test levels can skip registering it. Physics validation tests also get
  the `physics` label;
- functional scripts end with an assertion, not only a plot or a print.
