module.exports = {
  title: "MAGpurify2",
  base: "/magpurify2/",
  themeConfig: {
    nav: [
        { text: 'Repository', link: 'https://github.com/apcamargo/magpurify2' },
    ],
    displayAllHeaders: true,
    sidebar: {
      '/docs/': [
        '',
        'installation',
        'usage',
      ]
    }
  },
  plugins: [
    'vuepress-plugin-element-tabs',
  ],
}