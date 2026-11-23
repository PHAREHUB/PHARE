---
id: tests
name: Tests that prove the claim
paths:
  - "tests/**"
  - "pyphare/pyphare_tests/**"
---

For each test the PR adds or changes:
- would it fail if the bug came back or the feature were removed? If not, say
  which assertion is too weak;
- does it check the claim the PR makes, at the precision the PR claims? A
  tolerance much looser than the claimed precision is a finding;
- comparisons must fail on NaN (`EXPECT_NEAR` and `np.allclose` with NaN inputs,
  `assert x < tol` with x = NaN all pass silently in some forms);
- results compared bit for bit across runs are only reproducible with 1 or 2 MPI
  ranks (issue #1526); with 3 or more, compare with a tolerance;
- runtime: a test longer than about 10 minutes on the slowest PR build belongs in
  the nightly tier (`phare_test_labels(<id> nightly)` in CMake); physics
  validation tests also get the `physics` label;
- functional scripts end with an assertion, not only a plot or a print.
