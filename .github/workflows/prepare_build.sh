set -o errexit

case "$(uname -s)" in
    Linux)
        yum install -y clang zlib-devel bzip2-devel xz-devel
    ;;
    *)
esac

curl https://sh.rustup.rs -sSf | sh -s -- -y
rustup default nightly
pip install setuptools-rust