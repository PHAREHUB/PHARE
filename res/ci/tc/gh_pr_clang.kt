import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.XmlReport
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.buildFeatures.xmlReport
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhPrClang : BuildType({
    id("BuildGithubPrClang") // Phare_Phare_BuildGithubPrClang
    name = "gh_pr_clang"
    allowExternalStatus = true

    artifactRules = """
        **/gtest_out.xml => gtest_out.zip
        build/cppcheckHtml
    """.trimIndent()

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("install dependency", "gh_pr_clang", "00_install_dependency.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("configure", "gh_pr_clang", "01_configure.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("build", "gh_pr_clang", "02_build.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43")
        phareStep("run test", "gh_pr_clang", "03_run_test.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=5G")
        phareStep("run mpi test", "gh_pr_clang", "04_run_mpi_test.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=5G")
    }

    triggers {
        vcs {
            branchFilter = "+:*"
        }
    }

    failureConditions {
        executionTimeoutMin = 999
        errorMessage = true
    }

    features {
        xmlReport {
            reportType = XmlReport.XmlReportType.GOOGLE_TEST
            rules = "%system.teamcity.build.workingDir%/**/gtest_out.xml"
        }
        githubPullRequests()
        githubCommitStatus()
        perfmon {}
    }

    dependencies {
        artifacts(AbsoluteId("Phare_Phare_BuildSamraiUpstreamDep")) {
            id = "ARTIFACT_DEPENDENCY_4"
            enabled = false
            buildRule = lastSuccessful()
            cleanDestination = true
            artifactRules = "local.zip!** => %teamcity.build.workingDir%/.local-samrai"
        }
        artifacts(AbsoluteId("Phare_DistSamrai_Build")) {
            id = "ARTIFACT_DEPENDENCY_1"
            enabled = false
            buildRule = lastSuccessful()
            cleanDestination = false
            artifactRules = "local.zip!** => /root/.local-samrai"
        }
    }

    cleanup {
        keepLast(200)
    }

    requirements {
        equals("system.has_gcov", "true", "RQ_45")
        equals("system.has_lcov", "true", "RQ_46")
        equals("system.has_graphviz", "true", "RQ_49")
        equals("system.agent_name", "teamcity-docker-phare-fc31", "RQ_47")
        equals("env.TC_NODE_NAME", "phare-fc32-6cpu-20gb", "RQ_48")
        equals("env.N_CORES", "20")
    }

    disableSettings("RQ_45", "RQ_46", "RQ_49", "RQ_47", "RQ_48")
})
