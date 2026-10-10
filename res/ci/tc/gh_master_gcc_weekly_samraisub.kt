import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.triggers.ScheduleTrigger
import jetbrains.buildServer.configs.kotlin.triggers.schedule
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhMasterGccWeeklySamraisub : BuildType({
    id("GhMasterGccSamraisubNightly") // Phare_Phare_GhMasterGccSamraisubNightly
    name = "gh_master_gcc_weekly_samraisub"
    allowExternalStatus = true

    artifactRules = """
        data_out => data_out.zip
        **/gtest_out.xml => gtest_out.zip
        pyphare/pyphare_tests/*.png => pyphare_tests.zip
        build/tests/simulator/**/*.png => simulator_tests.zip
        build/tests/functional/**/*.png=>functional.png.zip
        build/tests/functional/**/*.mp4=>functional.mp4.zip
        build/tests/functional/**/*.pdf=>functional.pdf.zip
    """.trimIndent()

    params {
        password("sonarcloud_key", SONARCLOUD_KEY, display = ParameterDisplay.HIDDEN)
    }

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("config", "gh_master_gcc_weekly_samraisub", "00_config.sh", "129.104.6.165:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("build", "gh_master_gcc_weekly_samraisub", "01_build.sh", "129.104.6.165:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("test", "gh_master_gcc_weekly_samraisub", "02_test.sh", "129.104.6.165:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=12G")
        phareStep("testMPI_n_2", "gh_master_gcc_weekly_samraisub", "03_testmpi_n_2.sh", "129.104.6.165:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=12G")
        phareStep("caliper", "gh_master_gcc_weekly_samraisub", "04_caliper.sh", enabled = false)
    }

    triggers {
        vcs {
            enabled = false
            branchFilter = "+:*"
        }
        schedule {
            schedulingPolicy = weekly {
                dayOfWeek = ScheduleTrigger.DAY.Saturday
                hour = 1
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
        keepLast(200)
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
