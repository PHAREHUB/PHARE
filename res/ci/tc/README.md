# TeamCity pipelines

Step scripts for the PHARE pipelines on
[hephaistos TeamCity](https://hephaistos.lpp.polytechnique.fr/teamcity): `pr`, `nightly` and `weekly`.

## Layout

```
res/ci/tc/
  env.sh                 sourced by every step before its script
  <pipeline>/<job>/      one script per step, <NN>_<step>.sh
```

The pipelines are defined in TeamCity, each job step only does

```sh
. res/ci/tc/env.sh
. res/ci/tc/<pipeline>/<job>/<NN>_<step>.sh
```

so everything a step does is in its script here, and TeamCity only changes when jobs or steps are
added or removed. Scripts are spelled out per job rather than sharing helpers, so each one can be
read on its own. SAMRAI is always built as a subproject, so every job uses `teamcity-fedora`, and
tests run through [phlop](https://github.com/PhilipDeegan/phlop), serial and MPI together.

Steps of a job share the checkout, so a job's test step uses the `build/` of its build step.
Nothing is shared between jobs, which may run on different agents: every job builds what it tests.

## Pipelines

Docker run parameters below use `$N_CORES_ARG` for `-e N_CORES="$N_CORES"`.

### `pr`

Pull requests, all branches but the default one.

| job | image | step: docker run parameters |
|---|---|---|
| `gcc` | `.172` `teamcity-fedora:43` | `01_build`: `$N_CORES_ARG --shm-size=2G`, `02_test`: `$N_CORES_ARG --shm-size=5G`, `03_dry_run`: `$N_CORES_ARG --shm-size=1G` |
| `clang_asan` | `.172` `teamcity-fedora:43` | `01_build`: `$N_CORES_ARG --shm-size=2G`, `02_test`: `$N_CORES_ARG --shm-size=5G` |

`01_build` builds with `testMPI=OFF`, which builds every test target. `02_test` dumps the serial
test list, reconfigures with `-DtestMPI=ON` (only test registration changes) to dump the MPI one,
and runs both with phlop. `03_dry_run` sets up the functional tests without advancing them.

### `nightly`

Daily at 22:00 on the default branch.

| job | image | step: docker run parameters |
|---|---|---|
| `gcc` | `.165` `teamcity-fedora:43` | `01_build`: `$N_CORES_ARG --shm-size=2G`, `02_test`: `$N_CORES_ARG --shm-size=10G` |
| `clang_asan` | `.165` `teamcity-fedora:43` | `01_build`: `$N_CORES_ARG --shm-size=2G`, `02_test`: `$N_CORES_ARG --shm-size=10G` |

`01_build` builds two configurations, MPI tests up to exec level 100 (functional tests included) and
serial tests up to exec level 10, and dumps both test lists. `02_test` runs them together with phlop.
Artifacts, on the `02_test` step of both jobs:

```
.phlop => phlop.zip
data_out => data_out.zip
PHARE_REPORT.zip => PHARE_REPORT.zip
pyphare/pyphare_tests/**/*.png => pyphare_tests.zip
build/tests/simulator/**/*.png => simulator_tests.zip
build/tests/functional/**/*.png => functional.png.zip
build/tests/functional/**/*.mp4 => functional.mp4.zip
build/tests/functional/**/*.pdf => functional.pdf.zip
build/tests/functional/**/.phare => phare_stats.zip
build/tests/functional/**/.log => logs.zip
*.png => vtk_png.zip
```

### `weekly`

Sundays at 08:00 on the default branch.

| job | image | step: docker run parameters |
|---|---|---|
| `gcov` | `.165` `teamcity-fedora:44` | `01_build`: `$N_CORES_ARG --shm-size=1G`, `02_test`: `$N_CORES_ARG --shm-size=5G`, `03_gcovr`: `$N_CORES_ARG --shm-size=1G` |
| `caliper` | `.165` `teamcity-fedora:43` | `01_build`: `$N_CORES_ARG --shm-size=1G`, `02_test`: `$N_CORES_ARG --shm-size=12G` |

Artifacts: `gcov` `03_gcovr`: `build/coverage => coverage.zip`,
`caliper` `02_test`: `.phlop => phlop.zip`.

Release builds with serial and MPI tests, with and without ASan, run nightly, so the weekly
pipeline only has what nightly doesn't: a coverage build and caliper timings of the functional tests.

## Conventions

- Every step script starts with `set -ex`, if a script needs otherwise it should say why.
- Things every step needs (thread limits, `TMPDIR`, openmpi module, git `safe.directory`, etc) go in
  `env.sh`, things specific to a job stay in that job's scripts.
- TeamCity `%param%` references are not substituted inside these scripts, pass what's needed as
  `env.` parameters on the pipeline or job instead.

## Docker images

Every step runs in a container from the registry at `129.104.6.165:32219` or `129.104.6.172:32219`.

| image | SAMRAI | source repo | TeamCity build |
|---|---|---|---|
| `phare/teamcity-fedora:<fedora version>` | not installed, built as a subproject as needed | [phare-teamcity-agent](https://github.com/PHARCHIVE/phare-teamcity-agent) | [Build Docker images](https://hephaistos.lpp.polytechnique.fr/teamcity/admin/editBuild.html?id=buildType:Phare_BuildDockerImages) |
| `phare/teamcity-fedora_dep:<fedora version>` | installed into the container, use `-DSAMRAI_ROOT=/usr/local` (or `-DCMAKE_PREFIX_PATH`) | [phare-teamcity-agent_dep](https://github.com/PHARCHIVE/phare-teamcity-agent_dep) | [Build Docker-dep images](https://hephaistos.lpp.polytechnique.fr/teamcity/admin/editBuild.html?id=buildType:Phare_BuildDockerDepImages) |

The image builds run automatically when changes are merged into their source repo.
`Build Docker-dep images` also runs automatically when SAMRAI changes, so the SAMRAI
installed in `teamcity-fedora_dep` follows upstream.

`N_CORES` is set on each docker agent as the number of cores available to it,
and is passed into the containers via `-e N_CORES="$N_CORES"`, agents are selected
per job by the `env.N_CORES` requirement.
