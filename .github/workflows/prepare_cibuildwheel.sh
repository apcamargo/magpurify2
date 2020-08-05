# Abort on errors
set -e

# If on Linux, install dependencies for rust-htslib
case "$(uname -s)" in
    Linux)
        yum install -y clang zlib-devel bzip2-devel xz-devel
    ;;
    *)
esac

# Install Rust
curl https://sh.rustup.rs -sSf | sh -s -- -y --profile=minimal