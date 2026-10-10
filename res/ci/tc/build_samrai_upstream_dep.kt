import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.perfmon
import jetbrains.buildServer.configs.kotlin.buildSteps.ScriptBuildStep
import jetbrains.buildServer.configs.kotlin.triggers.vcs

object BuildSamraiUpstreamDep : BuildType({
    id("BuildSamraiUpstreamDep") // Phare_Phare_BuildSamraiUpstreamDep
    name = "build_samrai_upstream_dep"
    allowExternalStatus = true

    artifactRules = """
        %teamcity.build.workingDir%/samrai/build/docs/samrai-dox/html => doc.zip
    """.trimIndent()

    vcs {
        root(AbsoluteId("Phare_Phare_Samrai"))
        root(AbsoluteId(PHARE_VCS_ROOT), "+:res/ci/tc => .phare/res/ci/tc")
        cleanCheckout = true
    }

    steps {
        phareStep("load mpi,and build", "build_samrai_upstream_dep", "00_load_mpi_and_build.sh", "129.104.6.165:32219/phare/teamcity-fedora:43", ciRoot = ".phare/res/ci/tc", platform = ScriptBuildStep.ImagePlatform.Linux)
        phareStep("make doxygen doc", "build_samrai_upstream_dep", "01_make_doxygen_doc.sh", ciRoot = ".phare/res/ci/tc", enabled = false)
    }

    triggers {
        vcs {
            branchFilter = "+:*"
            triggerRules = "-:root=${PHARE_VCS_ROOT}:**"
        }
    }

    failureConditions {
        errorMessage = true
    }

    features {
        perfmon {}
    }

    dependencies {
        artifacts(AbsoluteId("Phare_Phare_BuildSamraiUpstreamDep")) {
            id = "ARTIFACT_DEPENDENCY_9"
            enabled = false
            buildRule = lastSuccessful()
            cleanDestination = false
            artifactRules = "ccache => %teamcity.build.workingDir%/.ccache"
        }
    }

    cleanup {
        keepLast(10)
    }

    requirements {
        equals("system.has_cppcheck", "true", "RQ_12")
        equals("system.has_gcov", "true", "RQ_13")
        equals("system.has_graphviz", "true", "RQ_14")
        equals("system.has_lcov", "true", "RQ_15")
        equals("env.N_CORES", "15")
    }

    disableSettings("RQ_12", "RQ_13", "RQ_14", "RQ_15")
})
