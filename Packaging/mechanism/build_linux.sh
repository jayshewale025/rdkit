#!/usr/bin/env bash
# Build, repair and verify the mechanism wheel inside a manylinux_2_28 container.
#   PY_TAG     CPython tag to build for (default: cp313), e.g. cp313, cp314
#   BUILD_DIR  build directory (default: build-mechanism-<arch>-<PY_TAG>)
#   JOBS       parallel compile jobs (default: number of CPUs)
# The architecture is the container's own (x86_64 or aarch64).
set -euo pipefail

arch="$(uname -m)"
case "$arch" in x86_64 | aarch64) ;; *) echo "unsupported architecture: $arch" >&2; exit 1 ;; esac
py_tag="${PY_TAG:-cp313}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build="${BUILD_DIR:-$root/build-mechanism-$arch-$py_tag}"
jobs="${JOBS:-$(nproc)}"
pyroot="/opt/python/$py_tag-$py_tag"
python="$pyroot/bin/python"
test -x "$python" || { echo "no CPython for $py_tag in this image ($python)" >&2; exit 1; }
export PATH="$pyroot/bin:$PATH"
echo "Building rdkit mechanism wheel: $py_tag / manylinux_2_28_$arch with $jobs jobs in $build"

# Fail before compilation if provenance cannot be read from the checkout.
git -C "$root" rev-parse --verify HEAD
git -C "$root" status --porcelain
dnf install -y cairo-devel freetype-devel bzip2-devel ccache
export CCACHE_DIR="$build/ccache"
ccache --max-size=2G
"$python" -m pip install -r "$root/Packaging/mechanism/requirements-build.txt" \
  cmake ninja 'numpy>=2,<3' Pillow pytest
mkdir -p "$build/dependencies" "$build/licenses"
cd "$build/dependencies"
curl --fail --location --retry 3 \
  https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.bz2 \
  --output boost_1_88_0.tar.bz2
echo "46d9d2c06637b219270877c9e16155cbd015b6dc84349af064c088e9b5b12f7b  boost_1_88_0.tar.bz2" | sha256sum --check -
tar -xf boost_1_88_0.tar.bz2
cd boost_1_88_0
./bootstrap.sh --with-python="$python" --with-libraries=python \
  --prefix="$build/boost"
if ! ./b2 -j"$jobs" variant=release link=shared threading=multi \
  cxxflags=-fPIC --with-python install > "$build/boost-build.log" 2>&1; then
  tail -n 120 "$build/boost-build.log"
  exit 1
fi
cp LICENSE_1_0.txt "$build/licenses/Boost-LICENSE_1_0.txt"
cp -a /usr/share/licenses "$build/licenses/system"
rpm -qa --qf '%{NAME}\t%{VERSION}-%{RELEASE}\t%{LICENSE}\t%{SOURCERPM}\n' \
  > "$build/licenses/system-packages.txt"
printf '%s\n' \
  'Boost source: https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.bz2' \
  'System source RPMs: https://vault.almalinux.org/ (AlmaLinux 8 repositories).' \
  'The system-packages inventory identifies the exact source RPM versions.' \
  > "$build/licenses/SOURCES.txt"

cmake -S "$root" -B "$build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_PREFIX_PATH="$build/boost" \
  -DPython_EXECUTABLE="$python" \
  -DRDK_BUILD_CPP_TESTS=OFF \
  -DBUILD_TESTING=OFF \
  -DRDK_BUILD_BOOST_PYTHON_WRAPPERS=ON \
  -DRDK_BUILD_NANOBIND_WRAPPERS=OFF \
  -DRDK_BUILD_CAIRO_SUPPORT=ON \
  -DRDK_BUILD_FREETYPE_SUPPORT=ON \
  -DRDK_INSTALL_COMIC_FONTS=OFF \
  -DRDK_BUILD_INCHI_SUPPORT=ON \
  -DRDK_INSTALL_INTREE=OFF \
  -DRDK_INSTALL_STATIC_LIBS=OFF \
  -DRDK_BUILD_COORDGEN_SUPPORT=OFF \
  -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
  -DRDK_BUILD_CHEMDRAW_SUPPORT=OFF \
  -DRDK_BUILD_PUBCHEMSHAPE_SUPPORT=OFF \
  -DRDK_BUILD_DESCRIPTORS3D=OFF \
  -DRDK_USE_BOOST_SERIALIZATION=OFF \
  -DRDK_USE_BOOST_IOSTREAMS=OFF \
  -DRDK_USE_BOOST_STACKTRACE=OFF
cmake --build "$build" --parallel "$jobs"
export LD_LIBRARY_PATH="$build/lib:$build/boost/lib"
"$python" "$root/Packaging/mechanism/build_wheel.py" \
  --build-dir "$build" --output "$build/wheelhouse" \
  --licenses-dir "$build/licenses" --platform "manylinux_2_28_$arch"

"$python" -m venv "$build/install-test"
venv="$build/install-test/bin/python"
"$venv" -m pip install "$build"/wheelhouse/*.whl pytest
mkdir -p "$build/installed-examples"
cd "$build/installed-examples"
unset LD_LIBRARY_PATH RDBASE PYTHONPATH
"$venv" -m rdkit_mechanism doctor
"$venv" -m rdkit_mechanism demo --output demo
"$venv" -m rdkit_mechanism gallery --png --output gallery
"$venv" -m rdkit_mechanism question6 --compact --output question6
"$venv" -m pytest -q --confcutdir="$root/Code/GraphMol/MolDraw2D/Wrap" \
  "$root/Code/GraphMol/MolDraw2D/Wrap/testMechanismAutomatic.py" \
  "$root/Code/GraphMol/MolDraw2D/Wrap/testMechanismDrawing.py" \
  "$root/Code/GraphMol/MolDraw2D/Wrap/testMechanismGallery.py"
