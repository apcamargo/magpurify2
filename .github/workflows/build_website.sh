# Abort on errors
set -e

# Install yarn and build website
yarn global add vuepress@1.5.3
yarn add -D maginapp/vuepress-plugin-katex
yarn --cwd ./website website:build
