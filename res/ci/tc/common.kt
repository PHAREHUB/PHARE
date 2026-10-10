import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.PullRequests
import jetbrains.buildServer.configs.kotlin.buildFeatures.commitStatusPublisher
import jetbrains.buildServer.configs.kotlin.buildFeatures.pullRequests
import jetbrains.buildServer.configs.kotlin.buildSteps.ScriptBuildStep
import jetbrains.buildServer.configs.kotlin.buildSteps.script

// One BuildType per <build>.kt in this directory, the step scripts live in <build>/<N>_<step>.sh
//  every step sources env.sh then its script from the checkout, so editing a step script
//  needs no TeamCity change, only adding/removing/reconfiguring steps does

const val PHARE_VCS_ROOT = "Phare_Phare_HttpsGithubComPharehubPhareRefsHeadsMaster"

// TODO: "credentialsJSON:<uuid>" references, created when Versioned Settings is enabled
//  the REST API / "View as code" do not expose the existing secure values
const val GITHUB_TOKEN = "credentialsJSON:TODO"
const val SONARCLOUD_KEY = "credentialsJSON:TODO"

const val N_CORES_ARG = "-e N_CORES=\"\$N_CORES\""

fun BuildSteps.phareStep(
    stepName: String,
    build: String,
    file: String,
    image: String? = null,
    runParams: String? = N_CORES_ARG,
    ciRoot: String = "res/ci/tc",
    platform: ScriptBuildStep.ImagePlatform? = null,
    enabled: Boolean = true,
    mode: BuildStep.ExecutionMode? = null,
) {
    script {
        name = stepName
        this.enabled = enabled
        mode?.let { executionMode = it }
        scriptContent = """
            . $ciRoot/env.sh
            . $ciRoot/$build/$file
        """.trimIndent()
        image?.let { dockerImage = it }
        platform?.let { dockerImagePlatform = it }
        if (image != null) runParams?.let { dockerRunParameters = it }
    }
}

fun BuildFeatures.githubPullRequests(enabled: Boolean = true) {
    pullRequests {
        this.enabled = enabled
        vcsRootExtId = PHARE_VCS_ROOT
        provider = github {
            authType = token { token = GITHUB_TOKEN }
            filterAuthorRole = PullRequests.GitHubRoleFilter.EVERYBODY
        }
    }
}

fun BuildFeatures.githubCommitStatus(enabled: Boolean = true) {
    commitStatusPublisher {
        this.enabled = enabled
        vcsRootExtId = PHARE_VCS_ROOT
        publisher = github {
            githubUrl = "https://api.github.com"
            authType = personalToken { token = GITHUB_TOKEN }
        }
    }
}

fun Cleanup.keepLast(count: Int, perBranch: Boolean = false) {
    keepRule {
        keepAtLeast = builds(count)
        dataToKeep = everything()
        applyPerEachBranch = perBranch
        preserveArtifactsDependencies = true
    }
}
