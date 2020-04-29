# Abort on errors
set -e

# Install yarn and build website
yarn global add vuepress@1.4.1
yarn --cwd ./website website:build
