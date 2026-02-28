# Releasing

## Branch Strategy

- **`dev`** — active development branch. Always carries a `-dev` pre-release version (e.g., `X.Y.Z-dev`).
- **`main`** — release-only branch. Each commit on `main` corresponds to a published version.
- Merge strategy: **squash merge** from `dev` into `main`.

## Release Checklist

1. **Update CHANGELOG.md**
   - Move items from `[Unreleased]` into a new `[X.Y.Z] - YYYY-MM-DD` section.
   - Add a comparison link at the bottom of the file.
   - Ensure all user-facing changes are documented.

2. **Finalize version**
   - Remove the `-dev` suffix from `version` in `Cargo.toml` (e.g., `X.Y.Z-dev` → `X.Y.Z`).

3. **Run checks**
   ```bash
   cargo check
   cargo fmt --check
   cargo clippy -- -D warnings
   cargo test
   ```

4. **Commit on `dev`**
   - Commit the CHANGELOG and version bump together.
   - Example: `Prepare release X.Y.Z`

5. **Squash merge into `main`**
   ```bash
   git checkout main
   git merge --squash dev
   git commit -m "Release X.Y.Z"
   ```

6. **Tag**
   ```bash
   git tag vX.Y.Z
   ```

7. **Push**
   ```bash
   git push origin main --tags
   ```

8. **Publish to crates.io**
   ```bash
   cargo publish --dry-run
   cargo publish
   ```

9. **Return to `dev` and bump to next dev version**
   ```bash
   git checkout dev
   git merge main
   ```
   - Bump `version` in `Cargo.toml` to the next anticipated version with `-dev` suffix.
   - Commit: `Start next development cycle (X.Y.Z-dev)`

## Versioning (SemVer)

This project follows [Semantic Versioning 2.0.0](https://semver.org/):

- **MAJOR** (`X.0.0`) — incompatible API changes.
- **MINOR** (`0.X.0`) — new functionality, backwards compatible.
- **PATCH** (`0.0.X`) — backwards-compatible bug fixes.

While the version is `0.x.y`, minor version bumps may include breaking changes.

## CHANGELOG Guidelines

Follow the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

### Categories

- **Added** — new features.
- **Changed** — changes to existing functionality.
- **Deprecated** — features that will be removed in upcoming releases.
- **Removed** — features that have been removed.
- **Fixed** — bug fixes.
- **Security** — vulnerability fixes.

### Rules

- Write entries from the user's perspective, not the developer's.
- Each entry should be a concise, complete sentence.
- Most recent release goes first.
- Always keep an `[Unreleased]` section at the top for ongoing work.
