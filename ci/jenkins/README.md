# Jenkins CI for Gkeyll

This sets up Jenkins (following the pattern used by
[SUNDIALS](https://github.com/llnl/sundials/tree/main/test/jenkins), adapted
for GitHub instead of Bitbucket) to build Gkeyll and run its unit tests on our
own persistent machines any time a pull request into `main` is opened or
updated.

Like `.github/workflows/mac_build.yml`, every build starts from a clean git
checkout and builds its dependencies (`gkylsoft/`) from scratch into its own
Jenkins workspace via `machines/mkdeps.macos.sh` — but on our own persistent
hardware instead of an ephemeral GitHub-hosted runner, and it actually
executes the unit tests rather than just building them.

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

## 5. Fix this node's `PATH`

Jenkins runs as a `launchd` service (via `brew services`), which does **not**
source your shell profile (`~/.zshrc`, `~/.zprofile`) — so it starts with
launchd's bare-bones default `PATH` (`/usr/bin:/bin:/usr/sbin:/sbin`) and
can't find Homebrew tools like `cmake`, even though `which cmake` works fine
in your terminal.

Fix this per-node in Jenkins itself (not by editing the launchd plist)
so it works the same way regardless of how each node's agent process is
started — the workstation you add later may be launched over SSH instead of
`brew services`, with a different environment altogether, so a plist edit on
this machine wouldn't help there anyway. Each node keeps its own value here,
which also means this generalizes cleanly once there's more than one node.

Manage Jenkins → Nodes → `manauref_lt1` → Configure → Node Properties → check
"Environment variables" → Add:

- Name: `PATH`
- Value: `/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin`

(Adjust if your Homebrew prefix or toolchain lives somewhere else.) When you
add the workstation later, give it the equivalent `PATH` for wherever *its*
build tools live.

## 6. Nothing to pre-build

Unlike an earlier version of this setup, the pipeline does **not** assume a
pre-built `gkylsoft/` on this machine. Each build passes
`PREFIX=$WORKSPACE/gkylsoft` to `machines/mkdeps.macos.sh` and
`machines/configure.macos.sh`, so dependencies are built from scratch inside
that build's own Jenkins workspace — isolated from any manual checkout you
have elsewhere, and safe for concurrent builds of different branches/PRs on
the same node (each gets its own workspace, so its own `gkylsoft/`).

The tradeoff is build time: expect each build to spend several extra minutes
building SuperLU/LuaJIT from source before it even gets to compiling Gkeyll,
same as `.github/workflows/mac_build.yml` does on GitHub's runners.

You still need Xcode command line tools / a working C toolchain and `cmake`
on this machine for `mkdeps.macos.sh` to succeed — if you've built Gkeyll
manually before, you already have these.

## 7. Create the multibranch pipeline job

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

## 8. Adding the second workstation later

1. On the workstation: nothing to install if using SSH launch — Jenkins pushes
   its agent jar over SSH. Make sure Java, a C/C++/Fortran toolchain, and
   `cmake` are present there — `mkdeps` builds everything else from scratch
   per build.
2. In Jenkins: Manage Jenkins → Nodes → New Node → "Permanent Agent" →
   Launch method "Launch agents via SSH", host/credentials for the
   workstation, label it e.g. `workstation-node`.
3. Same as step 5 above: add a `PATH` environment variable under that node's
   Node Properties covering wherever its build tools (`cmake`, compilers,
   MPI, etc.) actually live — don't assume it matches `manauref_lt1`'s value,
   especially if the workstation is Linux.
4. In `ci/jenkins/Jenkinsfile`, add an entry to the `nodes` map, using
   whichever `mkdeps`/`configure` scripts under `machines/` match that
   machine's OS (e.g. Linux instead of macOS):
   ```groovy
   def nodes = [
       'manauref_lt1': ['mkdeps.macos.sh', 'configure.macos.sh'],
       'workstation-node': ['mkdeps.linux.sh', 'configure.linux.cpu.sh'],
   ]
   ```
5. Commit and push. The next PR build will run on both machines in parallel.

## What the pipeline does

Per node, in parallel:

1. `printenv` — for debugging the build environment.
2. `make clean` and `rm -rf $WORKSPACE/gkylsoft` — the workspace persists
   across builds for a given branch/PR, so without this, a stale or (e.g.
   from a prior interrupted or racing build) corrupted file can look
   up-to-date and never get rebuilt/reinstalled, silently breaking the link
   or mixing old and new dependency files. Wiping both every run makes each
   build genuinely from scratch.
3. `machines/mkdeps.macos.sh` — builds `gkylsoft/` from scratch into
   `$WORKSPACE/gkylsoft`.
4. `machines/configure.macos.sh` — generates `config.mak` pointing at that
   freshly-built `gkylsoft/`.
5. `make -j3 check` — builds **and runs** all unit tests (`core`, `moments`,
   `vlasov`, `gyrokinetic`, `pkpm`); a failing unit test fails the build.
6. `make -j3 regression` — builds (does not execute) all regression tests,
   matching today's `.github/workflows/mac_build.yml` scope. There is no
   automated regression-test runner yet; running a curated subset of
   regression tests is a follow-up.
7. Archives `build/**/*.log` so logs are downloadable from the Jenkins build
   page even on failure.
