#!/usr/bin/env bash
# ====================================================================
# This file is part of PhaseTracer
#
# Installs the third-party dependencies of PhaseTracer.
#
#   System mode (default): installs packages with apt, dnf or brew (needs sudo on Linux).
#   Local mode (--local):  builds pinned sources (scripts/dependencies.cfg) into a prefix,
#                          default <repo>/.deps, as static -fPIC libraries (ALGLIB shared).
#                          No sudo needed; CMake picks the prefix up automatically.
# ====================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=dependencies.cfg
source "${SCRIPT_DIR}/dependencies.cfg"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Installs the dependencies of PhaseTracer: Boost (filesystem, log), ALGLIB, Eigen3, NLopt, GSL,
plus a C++ compiler, CMake and git.

Modes:
  (default)          Install system packages with apt, dnf or brew.
  --local            Build pinned sources into a prefix instead (no sudo). Boost, NLopt and GSL
                     are static, so PhaseTracer has no runtime dependency on them; ALGLIB is
                     shared (both PhaseTracer libraries use it) and found through the RPATH.

Options:
  --python           Also install what the Python interface needs (python3 headers, venv, pip).
  --hydrograv        Also pre-clone HydroGrav (branch ${HYDROGRAV_BRANCH}), so CMake needs no network.
  --prefix DIR       Prefix for --local (default: \$PT_DEPS_PREFIX or ${REPO_ROOT}/.deps).
  --only LIST        With --local, only build these (comma separated: ${LOCAL_DEPENDENCIES// /,}).
  --download-only    With --local, only download and verify the sources into <prefix>/src,
                     e.g. to copy them to a machine without internet access.
  --force            With --local, rebuild even if a dependency is already installed.
  -j, --jobs N       Parallel build jobs (default: number of CPUs).
  -y, --yes          Do not ask for confirmation (for CI).
  --dry-run          Print the commands instead of running them.
  -h, --help         Show this help.

After a --local install with a non-default prefix, configure with -DPT_DEPS_PREFIX=<prefix>.
EOF
}

# ------------------------------------------------------------------ helpers

log()  { printf '==> %s\n' "$*"; }
warn() { printf 'warning: %s\n' "$*" >&2; }
die()  { printf 'error: %s\n' "$*" >&2; exit 1; }

run() {
    printf '+ %s\n' "$*"
    if [[ ${DRY_RUN} -eq 0 ]]; then "$@"; fi
}

have() { command -v "$1" > /dev/null 2>&1; }

cpu_count() {
    if have nproc; then nproc
    elif have sysctl; then sysctl -n hw.ncpu
    else echo 2
    fi
}

sha256_of() {
    if have sha256sum; then sha256sum "$1" | cut -d' ' -f1
    elif have shasum; then shasum -a 256 "$1" | cut -d' ' -f1
    else die "need sha256sum or shasum to verify downloads"
    fi
}

# Upper-case variable prefix used in dependencies.cfg, e.g. nlopt -> NLOPT
cfg() { local var; var="$(echo "$1" | tr '[:lower:]' '[:upper:]')_$2"; echo "${!var}"; }

# ------------------------------------------------------------------ arguments

MODE=system
WITH_PYTHON=0
WITH_HG=0
DOWNLOAD_ONLY=0
FORCE=0
ASSUME_YES=0
DRY_RUN=0
ONLY=""
PREFIX="${PT_DEPS_PREFIX:-${REPO_ROOT}/.deps}"
JOBS="$(cpu_count)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --local)         MODE=local ;;
        --python)        WITH_PYTHON=1 ;;
        --hydrograv)     WITH_HG=1 ;;
        --prefix)        [[ $# -ge 2 ]] || die "--prefix needs a directory"; PREFIX="$2"; shift ;;
        --only)          [[ $# -ge 2 ]] || die "--only needs a list"; ONLY="$2"; shift ;;
        --download-only) DOWNLOAD_ONLY=1 ;;
        --force)         FORCE=1 ;;
        -j|--jobs)       [[ $# -ge 2 ]] || die "$1 needs a number"; JOBS="$2"; shift ;;
        -j*)             JOBS="${1#-j}" ;;
        -y|--yes)        ASSUME_YES=1 ;;
        --dry-run)       DRY_RUN=1 ;;
        -h|--help)       usage; exit 0 ;;
        *)               usage >&2; die "unknown option: $1" ;;
    esac
    shift
done

if [[ ${MODE} == system ]]; then
    [[ -z ${ONLY} && ${DOWNLOAD_ONLY} -eq 0 && ${FORCE} -eq 0 ]] \
        || die "--only, --download-only and --force only apply with --local"
fi

# ------------------------------------------------------------------ HydroGrav

clone_hydrograv() {
    local dest="${REPO_ROOT}/HydroGrav"
    if [[ -d ${dest} ]]; then
        log "HydroGrav already present in ${dest}"
        return
    fi
    have git || die "git is needed to clone HydroGrav"
    log "Cloning HydroGrav (${HYDROGRAV_BRANCH})"
    run git clone --branch "${HYDROGRAV_BRANCH}" --single-branch "${HYDROGRAV_URL}" "${dest}"
}

# ------------------------------------------------------------------ system mode

APT_PACKAGES="build-essential cmake git pkg-config libalglib-dev libnlopt-cxx-dev libeigen3-dev
              libboost-filesystem-dev libboost-log-dev libboost-thread-dev libgsl-dev"
APT_PYTHON="python3-dev python3-venv python3-pip"
DNF_PACKAGES="gcc-c++ make cmake git pkgconf alglib-devel nlopt-devel eigen3-devel boost-devel gsl-devel"
DNF_PYTHON="python3-devel python3-pip"
BREW_PACKAGES="cmake git pkg-config nlopt eigen boost gsl alglib libomp"
BREW_PYTHON="python"

install_system() {
    local sudo=""
    if [[ $(id -u) -ne 0 ]]; then
        if have sudo; then sudo="sudo"; else warn "not root and sudo not found; trying without"; fi
    fi
    local yes_flag=""
    [[ ${ASSUME_YES} -eq 1 ]] && yes_flag="-y"

    if have apt-get; then
        local pkgs="${APT_PACKAGES}"
        [[ ${WITH_PYTHON} -eq 1 ]] && pkgs="${pkgs} ${APT_PYTHON}"
        log "Installing packages with apt"
        # shellcheck disable=SC2086
        run ${sudo} apt-get update
        # shellcheck disable=SC2086
        run ${sudo} apt-get install ${yes_flag} ${pkgs}
    elif have dnf; then
        local pkgs="${DNF_PACKAGES}"
        [[ ${WITH_PYTHON} -eq 1 ]] && pkgs="${pkgs} ${DNF_PYTHON}"
        log "Installing packages with dnf"
        # shellcheck disable=SC2086
        run ${sudo} dnf install ${yes_flag} ${pkgs}
    elif have brew; then
        local pkgs="${BREW_PACKAGES}"
        [[ ${WITH_PYTHON} -eq 1 ]] && pkgs="${pkgs} ${BREW_PYTHON}"
        log "Installing packages with brew"
        # shellcheck disable=SC2086
        run brew install ${pkgs}
    else
        warn "no supported package manager found (apt-get, dnf, brew). Install the equivalents of:"
        printf '    %s\n' ${APT_PACKAGES} >&2
        die "or rerun with --local to build the libraries from source"
    fi
}

# ------------------------------------------------------------------ local mode

fetch() {
    local dep="$1"
    local file url sha dest
    file="$(cfg "${dep}" FILE)"; url="$(cfg "${dep}" URL)"; sha="$(cfg "${dep}" SHA256)"
    dest="${PREFIX}/src/${file}"

    if [[ -f ${dest} ]]; then
        [[ ${DRY_RUN} -eq 1 || $(sha256_of "${dest}") == "${sha}" ]] && return
        warn "checksum mismatch for cached ${file}; downloading again"
        rm -f "${dest}"
    fi

    log "Downloading ${dep} $(cfg "${dep}" VERSION)"
    if have curl; then
        run curl -fL --retry 3 -o "${dest}.part" "${url}"
    elif have wget; then
        run wget -O "${dest}.part" "${url}"
    else
        die "need curl or wget to download ${url}"
    fi
    [[ ${DRY_RUN} -eq 1 ]] && return

    local got
    got="$(sha256_of "${dest}.part")"
    if [[ ${got} != "${sha}" ]]; then
        rm -f "${dest}.part"
        die "checksum mismatch for ${file}: expected ${sha}, got ${got}"
    fi
    mv "${dest}.part" "${dest}"
}

# Extracts the archive of a dependency into a fresh build directory and echoes its path.
# Runs inside $(...), which does not inherit `set -e`, so failures are checked explicitly.
extract() {
    local dep="$1"
    local work="${PREFIX}/build"
    local dir="${work}/$(cfg "${dep}" DIR)"
    rm -rf "${dir}"
    tar --no-same-owner -xzf "${PREFIX}/src/$(cfg "${dep}" FILE)" -C "${work}" \
        || die "failed to extract $(cfg "${dep}" FILE)"
    [[ -d ${dir} ]] || die "expected ${dir} after extracting $(cfg "${dep}" FILE)"
    echo "${dir}"
}

cmake_install() {
    local src="$1"; shift
    run cmake -S "${src}" -B "${src}/_build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_PREFIX_PATH="${PREFIX}" \
        "$@"
    run cmake --build "${src}/_build" --parallel "${JOBS}"
    run cmake --install "${src}/_build"
}

build_eigen() {
    local src; src="$(extract eigen)" || exit 1
    cmake_install "${src}" -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_PKGCONFIG=ON
}

build_nlopt() {
    local src; src="$(extract nlopt)" || exit 1
    cmake_install "${src}" \
        -DBUILD_SHARED_LIBS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DNLOPT_PYTHON=OFF -DNLOPT_OCTAVE=OFF -DNLOPT_MATLAB=OFF -DNLOPT_GUILE=OFF \
        -DNLOPT_SWIG=OFF -DNLOPT_TESTS=OFF
}

build_gsl() {
    local src; src="$(extract gsl)" || exit 1
    (
        cd "${src}"
        run ./configure --prefix="${PREFIX}" --libdir="${PREFIX}/lib" \
            --disable-shared --enable-static --with-pic
        run make -j "${JOBS}"
        run make install
    )
}

# ALGLIB is linked by both libeffectivepotential and libphasetracer, so it is built shared:
# a static copy inside each library duplicates its global state and breaks it.
build_alglib() {
    local src; src="$(extract alglib)" || exit 1
    cp "${SCRIPT_DIR}/cmake/alglib/CMakeLists.txt" "${src}/CMakeLists.txt"
    rm -f "${PREFIX}"/lib/libalglib.*
    cmake_install "${src}" -DBUILD_SHARED_LIBS=ON
}

build_boost() {
    local src; src="$(extract boost)" || exit 1
    (
        cd "${src}"
        run ./bootstrap.sh --prefix="${PREFIX}" --with-libraries=log,filesystem,thread
        run ./b2 -j "${JOBS}" --prefix="${PREFIX}" --layout=system \
            link=static runtime-link=shared threading=multi variant=release \
            cxxflags=-fPIC cflags=-fPIC install
    )
}

check_local_tools() {
    local missing=""
    for tool in tar cmake make; do have "${tool}" || missing="${missing} ${tool}"; done
    have c++ || have g++ || have clang++ || missing="${missing} c++-compiler"
    have cc || have gcc || have clang || missing="${missing} c-compiler"
    have curl || have wget || missing="${missing} curl/wget"
    [[ -z ${missing} ]] || die "--local needs:${missing} (install them first, e.g. with the system mode)"
}

check_python_headers() {
    have python3 || { warn "python3 not found; the Python interface needs Python >= 3.8 with headers"; return; }
    if ! python3 -c 'import os, sysconfig; raise SystemExit(not os.path.exists(os.path.join(sysconfig.get_paths()["include"], "Python.h")))'; then
        warn "Python headers (Python.h) not found; install python3-dev / python3-devel from your system"
    fi
    python3 -c 'import venv' 2> /dev/null || warn "python3 venv module missing; install python3-venv"
}

check_openmp_macos() {
    [[ $(uname -s) == Darwin ]] || return 0
    if ! have brew || ! brew --prefix libomp > /dev/null 2>&1; then
        warn "Apple clang has no OpenMP runtime: run 'brew install libomp' (or the action loop runs serially)"
    fi
}

install_local() {
    local deps="${LOCAL_DEPENDENCIES}"
    if [[ -n ${ONLY} ]]; then
        deps="${ONLY//,/ }"
        for dep in ${deps}; do
            [[ " ${LOCAL_DEPENDENCIES} " == *" ${dep} "* ]] || die "unknown dependency '${dep}' (known: ${LOCAL_DEPENDENCIES})"
        done
    fi

    check_local_tools
    mkdir -p "${PREFIX}/src" "${PREFIX}/build"

    for dep in ${deps}; do fetch "${dep}"; done
    if [[ ${DOWNLOAD_ONLY} -eq 1 ]]; then
        log "Sources downloaded to ${PREFIX}/src; copy that directory and rerun with --local offline"
        return
    fi

    for dep in ${deps}; do
        local version stamp
        version="$(cfg "${dep}" VERSION)"
        stamp="${PREFIX}/${dep}.stamp"
        if [[ ${FORCE} -eq 0 && -f ${stamp} && $(cat "${stamp}") == "${version}" ]]; then
            log "${dep} ${version} already installed; skipping (use --force to rebuild)"
            continue
        fi
        if [[ ${DRY_RUN} -eq 1 ]]; then
            log "Would build ${dep} ${version} into ${PREFIX}"
            continue
        fi
        log "Building ${dep} ${version}"
        "build_${dep}"
        if [[ ${DRY_RUN} -eq 0 ]]; then
            echo "${version}" > "${stamp}"
            rm -rf "${PREFIX}/build/$(cfg "${dep}" DIR)"
        fi
    done

    [[ ${WITH_PYTHON} -eq 1 ]] && check_python_headers
    check_openmp_macos

    log "Dependencies installed in ${PREFIX}"
    if [[ ${PREFIX} != "${REPO_ROOT}/.deps" ]]; then
        log "Configure PhaseTracer with: cmake -DPT_DEPS_PREFIX=${PREFIX} ..."
    else
        log "CMake will pick them up automatically"
    fi
}

# ------------------------------------------------------------------ main

if [[ ${ASSUME_YES} -eq 0 && ${DRY_RUN} -eq 0 && -t 0 ]]; then
    if [[ ${MODE} == local ]]; then
        what="build dependencies into ${PREFIX}"
    else
        what="install system packages"
    fi
    read -r -p "This will ${what}. Continue? [y/N] " reply
    [[ ${reply} =~ ^[Yy]$ ]] || { log "Aborted"; exit 1; }
fi

if [[ ${MODE} == local ]]; then
    install_local
else
    install_system
fi

[[ ${WITH_HG} -eq 1 ]] && clone_hydrograv

log "Done"
