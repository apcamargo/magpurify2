# Abort on errors
set -e

# If on Linux, install dependencies for rust-htslib
case "$(uname -s)" in
    Linux)
        yum install -y clang zlib-devel bzip2-devel xz-devel
    ;;
    *)
esac

# Install Rust and set the toolchain to nightly
curl https://sh.rustup.rs -sSf | sh -s -- -y --profile=minimal
rustup default nightly

# Install setuptools-rust
pip install setuptools-rust==0.10.6