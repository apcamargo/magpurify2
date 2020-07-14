# Abort on errors
set -e

# Install yarn and build website
yarn global add vuepress@1.5.2
yarn add -D vuepress-plugin-mathjax
yarn --cwd ./website website:build
