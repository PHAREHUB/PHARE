# TeamCity builds

Build definitions and step scripts for the PHARE builds on
[hephaistos TeamCity](https://hephaistos.lpp.polytechnique.fr/teamcity) (`Phare / Phare` project).

## Layout

```
res/ci/tc/
  env.sh           sourced by every step before its script
  common.kt        shared Kotlin DSL helpers (phareStep, github features, keepLast)
  <build>.kt       one TeamCity BuildType per build, named as the build is on TeamCity
  <build>/         one script per build step, <NN>_<step name>.sh, NN = zero padded step order on TeamCity
```

Each TeamCity step only does

```sh
. res/ci/tc/env.sh
. res/ci/tc/<build>/<NN>_<step>.sh
```

from the checkout, so changing what a step does is just a change to its script here.
The `.kt` files only need changing when steps are added, removed, reordered or
reconfigured (docker image, run parameters, triggers, requirements, etc).

## Conventions

- Every step script starts with `set -ex`, if a script needs otherwise it should say why.
- Things every step needs (thread limits, `TMPDIR`, openmpi module, etc) go in `env.sh`,
  things specific to a build stay in that build's scripts.
- TeamCity `%param%` references are not substituted inside these scripts, pass what's
  needed as environment variables from the `.kt` instead.
- Disabled steps/settings are kept in the `.kt` files with `enabled = false` / `disableSettings`
  so the definitions match what is on TeamCity.

## Docker images

Every step runs in a container from the registry at `129.104.6.165:32219` or `129.104.6.172:32219`,
the image is set per step in the `.kt` files.

| image | SAMRAI | source repo | TeamCity build |
|---|---|---|---|
| `phare/teamcity-fedora:<fedora version>` | not installed, built as a subproject as needed | [phare-teamcity-agent](https://github.com/PHARCHIVE/phare-teamcity-agent) | [Build Docker images](https://hephaistos.lpp.polytechnique.fr/teamcity/admin/editBuild.html?id=buildType:Phare_BuildDockerImages) |
| `phare/teamcity-fedora_dep:<fedora version>` | installed into the container, use `-DSAMRAI_ROOT=/usr/local` (or `-DCMAKE_PREFIX_PATH`) | [phare-teamcity-agent_dep](https://github.com/PHARCHIVE/phare-teamcity-agent_dep) | [Build Docker-dep images](https://hephaistos.lpp.polytechnique.fr/teamcity/admin/editBuild.html?id=buildType:Phare_BuildDockerDepImages) |

The image builds run automatically when changes are merged into their source repo.
`Build Docker-dep images` also runs automatically when SAMRAI changes, so the SAMRAI
installed in `teamcity-fedora_dep` follows upstream.

Images currently used per build (enabled steps)

| build | image |
|---|---|
| `gh_pr_gcc`, `gh_pr_clang` | `.172` `teamcity-fedora_dep:43` |
| `gh_pr_gcc_samraisub`, `gh_pr_clang_samraisub` | `.172` `teamcity-fedora:43` |
| `gh_master_gcc_weekly_samraisub` | `.165` `teamcity-fedora:43` |
| `gh_master_gcc_samraisub_weekly_gcov` | `.165` `teamcity-fedora:44` |
| `gh_master_clang_weekly_asan` | `.165` `teamcity-fedora:42` |
| `gh_master_gcc_samraisub_weekly_caliper` | `.165` `teamcity-fedora_dep:43` |
| `gh_master_gcc_nightly_samraisub_phlop`, `gh_master_clang_nightly_samraisub_phlop_asan` | `.165` `teamcity-fedora_dep:43` |
| `build_samrai_upstream_dep` | `.165` `teamcity-fedora:43` |

`N_CORES` is set on each docker agent as the number of cores available to it,
and is passed into the containers via `-e N_CORES="$N_CORES"`, agents are selected
per build by the `env.N_CORES` requirement.

## Updating TeamCity

Not yet active: the `.kt` files are intended to be applied via TeamCity
Versioned Settings (Kotlin) on the `Phare_Phare` project. To get there

1. merge `res/ci/tc` to master, so branches being built have the step scripts
2. enable Versioned Settings on `Phare_Phare` with TeamCity committing its current settings
   (to a scratch branch), this generates the `credentialsJSON:<uuid>` references for the
   existing secure values
3. copy those references into `GITHUB_TOKEN` / `SONARCLOUD_KEY` in `common.kt`
4. point Versioned Settings at these files on master
