import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.triggers.ScheduleTrigger
import jetbrains.buildServer.configs.kotlin.triggers.schedule
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhMasterGccSamraisubWeeklyCaliper : BuildType({
    id("GhMasterGccSamraisubCaliperFunctional") // Phare_Phare_GhMasterGccSamraisubCaliperFunctional
    name = "gh_master_gcc_samraisub_weekly_caliper"
    allowExternalStatus = true

    artifactRules = """
        data_out => data_out.zip
        .phlop => phlop.zip
        **/gtest_out.xml => gtest_out.zip
        tests/simulator/*.png => simulator_tests.zip
        tests/functional/td/*.png => td.zip
        tests/functional/alfven_wave/*.png=>td.zip
    """.trimIndent()

    params {
        password("sonarcloud_key", SONARCLOUD_KEY, display = ParameterDisplay.HIDDEN)
    }

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("config", "gh_master_gcc_samraisub_weekly_caliper", "00_config.sh", "129.104.6.165:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("build", "gh_master_gcc_samraisub_weekly_caliper", "01_build.sh", "129.104.6.165:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("test", "gh_master_gcc_samraisub_weekly_caliper", "02_test.sh", "129.104.6.165:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=12G")
        phareStep("testMPI_n_2", "gh_master_gcc_samraisub_weekly_caliper", "03_testmpi_n_2.sh", enabled = false)
        phareStep("caliper", "gh_master_gcc_samraisub_weekly_caliper", "04_caliper.sh", enabled = false)
    }

    triggers {
        vcs {
            enabled = false
            branchFilter = "+:*"
        }
        schedule {
            schedulingPolicy = weekly {
                dayOfWeek = ScheduleTrigger.DAY.Sunday
                hour = 18
            }
            triggerBuild = always()
        }
    }

    failureConditions {
        executionTimeoutMin = 666
        errorMessage = true
    }

    features {
        githubPullRequests(enabled = false)
        githubCommitStatus(enabled = false)
        perfmon {}
    }

    cleanup {
        keepLast(55, perBranch = true)
    }

    requirements {
        equals("system.has_gcov", "true", "RQ_45")
        equals("system.has_lcov", "true", "RQ_46")
        equals("system.has_cppcheck", "true", "RQ_48")
        equals("system.has_graphviz", "true", "RQ_49")
        equals("system.agent_name", "teamcity-docker-phare-fc31", "RQ_47")
        equals("env.N_CORES", "44")
    }

    disableSettings("RQ_45", "RQ_46", "RQ_48", "RQ_49", "RQ_47")
})
