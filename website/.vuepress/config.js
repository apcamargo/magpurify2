module.exports = {
  title: "MAGpurify2",
  base: "/magpurify2/",
  themeConfig: {
    nav: [
        { text: 'Repository', link: 'https://github.com/apcamargo/magpurify2' },
    ],
    displayAllHeaders: true,
    sidebarDepth: 0,
    sidebar: {
      '/docs/': [
        '',
        {
          title: 'Theory',
          collapsable: false,
          children: [
            'composition-theory',
            'coverage-theory',
            'codon-usage-theory',
          ]
        },
        'installation',
        'usage',
      ]
    }
  },
  plugins: [
    [
      '@maginapp/vuepress-plugin-katex',
        {
          delimiters: 'dollars'
        }
    ],
  ],
}