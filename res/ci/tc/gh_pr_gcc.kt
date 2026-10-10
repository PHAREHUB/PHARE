import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.XmlReport
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.buildFeatures.xmlReport
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object GhPrGcc : BuildType({
    id("BuildGithubPr") // Phare_Phare_BuildGithubPr
    name = "gh_pr_gcc"
    allowExternalStatus = true
    publishArtifacts = PublishMode.ALWAYS

    artifactRules = """
        **/gtest_out.xml => gtest_out.zip
        coverage => coverage.zip
        subprojects/pharead/doc => documentation.zip
        cppcheckHtml => cppcheck.zip
        stats => gitstats.zip
        html => dev-doc.zip
        tests/simulator/*.png => simulator_tests.zip
    """.trimIndent()

    params {
        password("sonarcloud_key", SONARCLOUD_KEY, display = ParameterDisplay.HIDDEN)
    }

    vcs {
        root(AbsoluteId(PHARE_VCS_ROOT))
        cleanCheckout = true
    }

    steps {
        phareStep("install dependency", "gh_pr_gcc", "00_install_dependency.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=1G")
        phareStep("configure", "gh_pr_gcc", "01_configure.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=2G")
        phareStep("Create SonarQube properties file", "gh_pr_gcc", "02_create_sonarqube_properties_file.sh", enabled = false)
        phareStep("build", "gh_pr_gcc", "03_build.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43")
        phareStep("pre_coverage", "gh_pr_gcc", "04_pre_coverage.sh", enabled = false)
        phareStep("run test", "gh_pr_gcc", "05_run_test.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=5G")
        phareStep("run MPI test", "gh_pr_gcc", "06_run_mpi_test.sh", "129.104.6.172:32219/phare/teamcity-fedora_dep:43", "$N_CORES_ARG --shm-size=1G")
        step {
            name = "Vote on Pull Request"
            type = "VoteRhodecodePr"
            enabled = false
            executionMode = BuildStep.ExecutionMode.RUN_ON_FAILURE
            param("PR_FAILED_MESSAGE", """
                Teamcity bot:

                <a href="%teamcity.serverUrl%/viewLog.html?buildId=%teamcity.build.id%">
                <img src="%teamcity.serverUrl%/app/rest/builds/id:%teamcity.build.id%/statusIcon.svg?guest=1">
                </a>

                Your PR failed in one ore more build steps.

                [Build overview](%teamcity.serverUrl%/viewLog.html?buildId=%teamcity.build.id%&tab=buildResultsDiv&guest=1)
                [Build logs](%teamcity.serverUrl%/viewLog.html?buildId=%teamcity.build.id%&tab=buildLog&guest=1)
            """.trimIndent())
            param("PR_SUCCESS_MESSAGE", """
                Teamcity bot:

                <a href="%teamcity.serverUrl%/viewLog.html?buildId=%teamcity.build.id%">
                <img src="%teamcity.serverUrl%/app/rest/builds/id:%teamcity.build.id%/statusIcon.svg?guest=1">
                </a>

                Congratulations, your PR passes all tests!

                [Build overview](%teamcity.serverUrl%/viewLog.html?buildId=%teamcity.build.id%&tab=buildResultsDiv&guest=1)
                [Build logs](%teamcity.serverUrl%/viewLog.html?buildId=%teamcity.build.id%&tab=buildLog&guest=1)
            """.trimIndent())
        }
        phareStep("coverage_report", "gh_pr_gcc", "07_coverage_report.sh", enabled = false)
        phareStep("coverage_report (1)", "gh_pr_gcc", "08_coverage_report_1.sh", enabled = false)
        phareStep("documentation", "gh_pr_gcc", "09_documentation.sh", enabled = false, mode = BuildStep.ExecutionMode.RUN_ON_FAILURE)
        phareStep("cppcheck", "gh_pr_gcc", "10_cppcheck.sh", enabled = false)
        phareStep("generate gitstats", "gh_pr_gcc", "11_generate_gitstats.sh", "129.104.6.165:32219/phare/teamcity-fedora_dep:37", "$N_CORES_ARG --shm-size=1G", enabled = false, mode = BuildStep.ExecutionMode.RUN_ON_FAILURE)
        phareStep("Publish on SonarCloud", "gh_pr_gcc", "12_publish_on_sonarcloud.sh", enabled = false, mode = BuildStep.ExecutionMode.RUN_ON_FAILURE)
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
            id = "ARTIFACT_DEPENDENCY_3"
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
        equals("system.has_cppcheck", "true", "RQ_48")
        equals("system.has_graphviz", "true", "RQ_49")
        equals("system.agent_name", "teamcity-docker-phare-fc31", "RQ_47")
        equals("env.N_CORES", "30")
    }

    disableSettings("RQ_45", "RQ_46", "RQ_48", "RQ_49", "RQ_47")
})
