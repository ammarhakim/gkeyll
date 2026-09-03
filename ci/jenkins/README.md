# Jenkins CI for Gkeyll

This sets up Jenkins (following the pattern used by
[SUNDIALS](https://github.com/llnl/sundials/tree/main/test/jenkins), adapted
for GitHub instead of Bitbucket) to build Gkeyll and run its unit tests on our
own persistent machines any time a pull request into `main` is opened or
updated.

Unlike `.github/workflows/mac_build.yml`, which spins up a fresh, empty
GitHub-hosted runner for every build, these Jenkins nodes are long-lived
machines that already have the `gkylsoft/` toolchain installed. The pipeline
therefore skips dependency installation (`install-deps/mkdeps.sh`) and goes
straight to building and testing.

The controller runs on `manauref_lt1` (this laptop). A second workstation can be
added later as an additional build node — see the last section.

## 1. Install Jenkins

```
brew install jenkins-lts
brew services start jenkins-lts
```

This pulls in `openjdk@21` as a dependency. Jenkins listens on
`http://localhost:8080`. Get the initial admin password with:

```
cat ~/.jenkins/secrets/initialAdminPassword
```

(`~/.jenkins` is Homebrew's default `JENKINS_HOME` on macOS — not
`/opt/homebrew/var/lib/jenkins`, despite what some older docs say.) Open
`http://localhost:8080`, paste the password, and install the "suggested
plugins" set when prompted.

## 2. Install additional plugins

Manage Jenkins → Plugins → Available plugins, install:

- **GitHub Branch Source** — lets Jenkins discover branches/PRs on a GitHub
  repo and is the GitHub equivalent of the Bitbucket Branch Source plugin
  SUNDIALS uses.
- **Collapsing Console Sections** (optional) — collapses long build log
  sections for readability, same as SUNDIALS' setup.

## 3. Add a GitHub credential

Jenkins needs read access to `gkeyllorg/gkeyll` to poll for branches and PRs:

1. Create a GitHub Personal Access Token with `repo` scope (Settings →
   Developer settings → Personal access tokens on GitHub).
2. In Jenkins: Manage Jenkins → Credentials → System → Global credentials →
   Add Credentials → kind "GitHub personal access token" (or "Username with
   password", username = your GitHub username, password = the token).

## 4. Label this node `manauref_lt1`

Manage Jenkins → Nodes → "Built-In Node" → Configure → Labels: add
`manauref_lt1`. This is the label the Jenkinsfile's `node('manauref_lt1')` matches.

## 5. Confirm the toolchain is already built

The pipeline assumes `gkylsoft/` already exists on this machine via the usual
flow (from the `gkeyll/` directory):

```
../gkylsoft  # should already exist
./machines/mkdeps.macos.sh   # one-time, only if gkylsoft/ is missing
./machines/configure.macos.sh
make -j3 install
```

If this hasn't been done yet, do it once manually before relying on Jenkins.

## 6. Create the multibranch pipeline job

New Item → name it (e.g. `gkeyll`) → **Multibranch Pipeline**.

- **Branch Sources** → Add source → GitHub:
  - Credentials: the one added in step 3.
  - Repository HTTPS URL: `https://github.com/gkeyllorg/gkeyll`
  - Behaviors: "Discover branches" (at least `main`) and "Discover pull
    requests from origin". Add "Discover pull requests from forks" only if
    external contributors submit PRs from forks.
- **Build Configuration** → Mode: "by Jenkinsfile", Script Path:
  `ci/jenkins/Jenkinsfile`.
- **Scan Multibranch Pipeline Triggers** → check "Periodically if not
  otherwise run" → interval **2 minutes**.

Since this laptop is behind NAT, GitHub can't push a webhook to it. The
periodic scan above is how Jenkins learns about new/updated PRs instead —
every 2 minutes it asks GitHub for changes and kicks off a build if there are
any. This means there's up to a ~2 minute delay before a build starts, which
is an acceptable tradeoff for not having to run a public tunnel (e.g. ngrok)
just to receive webhooks.

Save. Jenkins will scan the repo, find `main` and any open PRs with a
Jenkinsfile, and start building them.

## 7. Adding the second workstation later

1. On the workstation: nothing to install if using SSH launch — Jenkins pushes
   its agent jar over SSH. Make sure Java and the `gkylsoft/` toolchain are
   present there (steps 5, mirrored on that machine).
2. In Jenkins: Manage Jenkins → Nodes → New Node → "Permanent Agent" →
   Launch method "Launch agents via SSH", host/credentials for the
   workstation, label it e.g. `workstation-node`.
3. In `ci/jenkins/Jenkinsfile`, uncomment/add a line to the `nodeLabels` map:
   ```groovy
   def nodeLabels = [
       'manauref_lt1': 'manauref_lt1',
       'workstation-node': 'workstation-node',
   ]
   ```
4. Commit and push. The next PR build will run on both machines in parallel.

## What the pipeline does

Per node, in parallel:

1. `printenv` — for debugging the build environment.
2. `make -j3 check` — builds **and runs** all unit tests (`core`, `moments`,
   `vlasov`, `gyrokinetic`, `pkpm`); a failing unit test fails the build.
3. `make -j3 regression` — builds (does not execute) all regression tests,
   matching today's `.github/workflows/mac_build.yml` scope. There is no
   automated regression-test runner yet; running a curated subset of
   regression tests is a follow-up.
4. Archives `build/**/*.log` so logs are downloadable from the Jenkins build
   page even on failure.
