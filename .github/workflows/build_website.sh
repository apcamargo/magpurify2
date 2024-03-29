# Abort on errors
set -e

# Install yarn and build website
yarn global add vuepress@1.8.2
yarn add -D maginapp/vuepress-plugin-katex markdown-it-footnote
yarn --cwd ./website website:build
