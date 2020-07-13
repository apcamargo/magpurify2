# Abort on errors
set -e

# Install yarn and build website
yarn add -D vuepress@1.5.2 vuepress-plugin-md-enhance@0.7.0
yarn --cwd ./website website:build
