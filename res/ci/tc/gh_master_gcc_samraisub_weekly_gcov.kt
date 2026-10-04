import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.triggers.ScheduleTrigger
import jetbrains.buildServer.configs.kotlin.triggers.schedule
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhMasterGccSamraisubWeeklyGcov : BuildType({
    id("GhMasterGccSamraisubNightlyGcov") // Phare_Phare_GhMasterGccSamraisubNightlyGcov
    name = "gh_master_gcc_samraisub_weekly_gcov"
    description = "Weekly on Sunday at 08:00"
    allowExternalStatus = true

    artifactRules = """
        **/gtest_out.xml => gtest_out.zip
        build/coverage => coverage.zip
    """.trimIndent()

    params {
        password("sonarcloud_key", SONARCLOUD_KEY, display = ParameterDisplay.HIDDEN)
    }

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("config", "gh_master_gcc_samraisub_weekly_gcov", "00_config.sh", "129.104.6.165:32219/phare/teamcity-fedora:44")
        phareStep("build", "gh_master_gcc_samraisub_weekly_gcov", "01_build.sh", "129.104.6.165:32219/phare/teamcity-fedora:44")
        phareStep("test", "gh_master_gcc_samraisub_weekly_gcov", "02_test.sh", "129.104.6.165:32219/phare/teamcity-fedora:44", "$N_CORES_ARG --shm-size=5G")
        phareStep("gcovr", "gh_master_gcc_samraisub_weekly_gcov", "03_gcovr.sh", "129.104.6.165:32219/phare/teamcity-fedora:44")
    }

    triggers {
        vcs {
            enabled = false
            branchFilter = "+:*"
        }
        schedule {
            schedulingPolicy = weekly {
                dayOfWeek = ScheduleTrigger.DAY.Sunday
                hour = 8
            }
            triggerBuild = always()
        }
    }

    failureConditions {
        executionTimeoutMin = 2222
        errorMessage = true
    }

    features {
        githubPullRequests(enabled = false)
        githubCommitStatus(enabled = false)
        perfmon {}
    }

    cleanup {
        keepLast(50)
    }

    requirements {
        equals("system.has_gcov", "true", "RQ_45")
        equals("system.has_lcov", "true", "RQ_46")
        equals("system.has_cppcheck", "true", "RQ_48")
        equals("system.has_graphviz", "true", "RQ_49")
        equals("system.agent_name", "teamcity-docker-phare-fc31", "RQ_47")
        equals("env.N_CORES", "20")
    }

    disableSettings("RQ_45", "RQ_46", "RQ_48", "RQ_49", "RQ_47")
})
