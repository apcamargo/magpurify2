# abort on errors
set -e

# install yarn and build
yarn global add vuepress@1.4.1
yarn --cwd ./website website:build
