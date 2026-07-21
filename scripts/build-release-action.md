# Meryl release with GitHub Actions

## 1. Modify versions

* ../src/main.mk
  ```
  VERSION      := release v1.4.2 # vX.X.X, needs `v`
  ```

* version_update.pl
  ```
  my $label    = "release";      #  If not 'release' print this in the version output. Default: 'snapshot'.
  my $major    = "1";            #  Bump before release.
  my $minor    = "4.2";          #  Bump before release.
  ```

* README.md
  * Add Change log
  * Bump version-specific links

## 2. Commit and tag
* Trigger Actions workflow with a tag
  ```sh
  git add src/main.mk scripts/version_update.pl README.md
  git commit -m "release v1.4.2"
  git push
  git tag -a v1.4.2 -m "meryl 1.4.2"
  git push origin v1.4.2
  ```

* In case we need to revert, remove tag
  ```sh
  git tag -d v1.4.2
  git push origin :refs/tags/v1.4.2
  ```

## 3. Test built binaries
  * `./bin/meryl --version`
  * Do some test runs with small input files

## 4. Publish
  * Actions generates a draft release with Assets
  * Update release notes
  * Publish

## 5. Restore back to snapshot

* ../src/main.mk
  ```
  VERSION      := snapshot 1.4.2 # X.X.X
  ```

* version_update.pl
  ```
  my $label    = "snapshot";     #  If not 'release' print this in the version output. Default: 'snapshot'.
  my $major    = "1";            #  Bump before release.
  my $minor    = "4.2";          #  Bump before release.
  ```

  ```sh
  git commit -m "post-release: bump to snapshot 1.4.2"
  ```

## MacOS notes

The release workflow builds two macOS tarballs:

- `Darwin-arm64` on `macos-14` (Apple Silicon)
- `Darwin-amd64` on `macos-15-intel` (Intel, macOS 15 — GitHub commits
  to supporting this runner through Fall 2027)
- See [notes](https://github.blog/changelog/2025-09-19-github-actions-macos-13-runner-image-is-closing-down/) from GitHub

Both tarballs are ad-hoc signed. That is sufficient for a personal Mac
where the user runs `xattr -dr com.apple.quarantine` on the extracted
tree, but it is **not** sufficient on MDM-managed Macs — Gatekeeper
will kill the binaries regardless of architecture. Users on locked-down
machines should build from source (see below); a proper fix requires
Apple Developer ID signing and notarization, which is out of scope here.

## Building meryl locally

Meryl builds with plain `gmake` on Linux and macOS. No CMake, no
autotools.

### Prerequisites

- `git` 2.12+ for submodules
- GNU make 3.81+ (`gmake` on macOS; `make` on Linux)
- A C++17 compiler (gcc 9+ or clang 12+)
- Dev headers for: `zlib`, `bzip2`, `xz` (liblzma), `openssl`, `libcurl`

Install prerequisites:

- **Ubuntu / Debian**
  ```sh
  sudo apt install build-essential git \
      zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev
  ```
- **RHEL / Rocky / Alma / Fedora**
  ```sh
  sudo dnf install gcc-c++ make git \
      zlib-devel bzip2-devel xz-devel openssl-devel libcurl-devel
  ```
- **macOS** (Xcode Command Line Tools + Homebrew)
  ```sh
  xcode-select --install
  brew install make gcc openssl libcurl xz
  ```

### Build

```sh
git clone --recursive https://github.com/marbl/meryl
cd meryl/src
gmake -j      # on macOS: gmake, on Linux: make (both are GNU make)
```

Binaries are placed in `../build/bin/`:

```sh
../build/bin/meryl --version
../build/bin/meryl --help
```

To install to a system path, copy `../build/bin/*` and `../build/lib/*`
somewhere on your `PATH`, or set `PATH="$PWD/../build/bin:$PATH"`.

### Notes

- The `--recursive` clone flag pulls in the `meryl-utility` submodule
  under `src/utility`. Without it, the build fails with
  `Failed to create version.H`. If you already cloned without
  `--recursive`, run `git submodule update --init --recursive`.
- Rebuild in place after edits with `gmake -j`; use `gmake clean` to
  wipe `build/`.
- macOS binaries built locally are ad-hoc signed by the linker and run
  fine on your own machine. They will NOT satisfy Gatekeeper if
  redistributed via download — that requires Apple Developer ID
  signing and notarization, which is out of scope here.
