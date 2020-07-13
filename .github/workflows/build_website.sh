# Abort on errors
set -e

# Install yarn and build website
yarn global add vuepress@1.5.2
yarn add -D vuepress-plugin-md-enhance@0.7.0
yarn --cwd ./website website:build
