# Abort on errors
set -e

# Install yarn and build website
yarn global add vuepress@1.6.0
yarn add -D maginapp/vuepress-plugin-katex markdown-it-footnote
yarn --cwd ./website website:build
