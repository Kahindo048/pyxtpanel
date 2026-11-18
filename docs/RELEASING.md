# Releasing pyxtpanel

A release is made by pushing a tag. GitHub runs the full test matrix, builds the package, publishes it on
PyPI and creates the GitHub release; nothing is uploaded from a personal computer, and no password or
token is stored anywhere (PyPI "trusted publishing"). The workflows are in `.github/workflows/`:
`ci.yml` (checks on every push and pull request, called again by the release) and `release.yml`.

## One-time setup

1. **The GitHub repository.** Create `pyxtpanel` under the account that owns the other packages and push
   `main`. The links of the README and `pyproject.toml` assume `github.com/Kahindo048/pyxtpanel`.
2. **The environments.** In the repository: *Settings, Environments*, create `pypi` and `testpypi`. On
   `pypi`, add yourself as a *required reviewer*: a tag then waits for one click before anything is
   published.
3. **Trusted publisher on PyPI.** The project does not exist yet, so use a *pending publisher*: on
   pypi.org, *Your account, Publishing, Add a new pending publisher*, with project name `pyxtpanel`,
   owner, repository `pyxtpanel`, workflow `release.yml`, environment `pypi`. Do the same on
   test.pypi.org (a separate account) with environment `testpypi`.
4. **Protect `main`.** *Settings, Branches*: require the `CI` checks (`Quality`, `Tests`, `Build`) to
   pass before merging.

## Each release

1. **Version.** Change it in the three places that must agree (`tests/contract/test_release_files.py`
   fails if they differ): `pyproject.toml`, `CITATION.cff`, `src/pyxtpanel/_version.py`.
2. **Changelog.** In `CHANGELOG.md`, rename `[Unreleased]` to `[X.Y.Z] - date` and open a new
   `[Unreleased]`.
3. **Locally**, as in CONTRIBUTING: `ruff check . && ruff format --check . && mypy && lint-imports`,
   `pytest`, `python -m tools.doc --write`, then `python -m build` and `python -m twine check dist/*`.
4. **Rehearsal on TestPyPI.** *Actions, Release, Run workflow* (repository `testpypi`). Then, in a fresh
   environment: `pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/
   pyxtpanel==X.Y.Z` and run `python -I tools/smoke_install.py`. A version can be uploaded to an index
   only once, so a rehearsal uses up the number on TestPyPI: use a pre-release such as `0.1.0rc1`.
5. **Tag.** `git tag vX.Y.Z` then `git push origin vX.Y.Z`. The tag must be the version of
   `pyproject.toml` with a `v` in front; the workflow refuses it otherwise. A tag with `rc`, `a`, `b` or
   `dev` in it creates a GitHub *pre-release*.
6. **Check.** The PyPI page (the README is its description), `pip install pyxtpanel` in a clean
   environment, and the GitHub release with the two files attached.

## If something goes wrong

* **A red build.** Nothing is published: `release.yml` starts with the whole CI.
* **A bad release already on PyPI.** A file cannot be replaced or re-uploaded under the same version.
  *Yank* the release on PyPI (it stays installable by exact version, but is no longer chosen by
  `pip install pyxtpanel`) and publish the next number.
* **The publish job is refused** (`invalid-publisher`): the workflow name, the environment or the
  repository in the trusted-publisher form on PyPI does not match `release.yml`, `pypi` and the
  repository that ran it.

## The estimation packages

`pyxtabond2`, `pyxthreg` and `pyxtdpthresh` can copy `ci.yml` and `release.yml` as they are (change the
project name in the two `url:` lines and the paths of the package). Once they depend on this package,
pin it as `pyxtpanel>=0.3,<2` (see the contract in part 8 of the documentation).
