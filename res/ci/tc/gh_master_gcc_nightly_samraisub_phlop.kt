import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.triggers.schedule
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhMasterGccNightlySamraisubPhlop : BuildType({
    id("GhMasterGccSamraisubNightlyPhlop") // Phare_Phare_GhMasterGccSamraisubNightlyPhlop
    name = "gh_master_gcc_nightly_samraisub_phlop"
    allowExternalStatus = true

    artifactRules = """
        PHARE_REPORT.zip => PHARE_REPORT.zip
        data_out => data_out.zip
        .phlop => phlop.zip
        pypyphare_tests/**/*.png => pyphare_tests.zip
        build/tests/simulator/**/*.png => simulator_tests.zip
        build/tests/functional/**/*.png=>functional.png.zip
        build/tests/functional/**/*.mp4=>functional.mp4.zip
        build/tests/functional/**/*.pdf=>functional.pdf.zip
        build/tests/functional/**/.phare=>phare_stats.zip
        build/tests/functional/**/.log=>logs.zip
        *.png=>vtk_png.zip
    """.trimIndent()

    params {
        password("sonarcloud_key", SONARCLOUD_KEY, display = ParameterDisplay.HIDDEN)
    }

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("run", "gh_master_gcc_nightly_samraisub_phlop", "00_run.sh", "129.104.6.165:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=10G")
    }

    triggers {
        vcs {
            enabled = false
            branchFilter = "+:*"
        }
        schedule {
            schedulingPolicy = daily { hour = 22 }
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
