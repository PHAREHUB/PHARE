import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.XmlReport
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.buildFeatures.xmlReport
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhPrClangSamraisub : BuildType({
    id("BuildGithubPrClangSs") // Phare_Phare_BuildGithubPrClangSs
    name = "gh_pr_clang_samraisub"
    description = "build with SAMRAI as subproject"
    allowExternalStatus = true

    artifactRules = """
        **/gtest_out.xml => gtest_out.zip
        build/tests/simulator/**/*.png => simulator_tests.zip
    """.trimIndent()

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("install dependency", "gh_pr_clang_samraisub", "00_install_dependency.sh", "129.104.6.172:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("configure", "gh_pr_clang_samraisub", "01_configure.sh", "129.104.6.172:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("build", "gh_pr_clang_samraisub", "02_build.sh", "129.104.6.172:32219/phare/teamcity-fedora:43")
        phareStep("test", "gh_pr_clang_samraisub", "03_test.sh", "129.104.6.172:32219/phare/teamcity-fedora:43", "$N_CORES_ARG --shm-size=1G")
    }

    triggers {
        vcs {
            branchFilter = "+:*"
        }
    }

    failureConditions {
        executionTimeoutMin = 333
        errorMessage = true
    }

    features {
        xmlReport {
            enabled = false
            reportType = XmlReport.XmlReportType.GOOGLE_TEST
            rules = "%system.teamcity.build.workingDir%/**/gtest_out.xml"
        }
        githubPullRequests()
        githubCommitStatus()
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
        equals("env.N_CORES", "20")
    }

    disableSettings("RQ_45", "RQ_46", "RQ_48", "RQ_49", "RQ_47")
})
