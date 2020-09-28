(window.webpackJsonp=window.webpackJsonp||[]).push([[13],{372:function(a,t,s){"use strict";s.r(t);var e=s(17),r=Object(e.a)({},(function(){var a=this,t=a.$createElement,s=a._self._c||t;return s("ContentSlotsDistributor",{attrs:{"slot-key":a.$parent.slotKey}},[s("h1",{attrs:{id:"installation"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#installation"}},[a._v("#")]),a._v(" Installation")]),a._v(" "),s("h2",{attrs:{id:"installing-magpurify2"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#installing-magpurify2"}},[a._v("#")]),a._v(" Installing MAGpurify2")]),a._v(" "),s("p",[a._v("MAGpurify2 can be installed using "),s("code",[a._v("pip")]),a._v(" or "),s("code",[a._v("conda")]),a._v(". Alternatively, you can execute it via "),s("RouterLink",{attrs:{to:"/docs/installation/#docker"}},[a._v("Docker")]),a._v(".")],1),a._v(" "),s("code-group",[s("code-block",{attrs:{title:"pip",active:""}},[s("div",{staticClass:"language-bash extra-class"},[s("pre",{pre:!0,attrs:{class:"language-bash"}},[s("code",[a._v("pip "),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("install")]),a._v(" magpurify2\n")])])])]),a._v(" "),s("code-block",{attrs:{title:"conda"}},[s("div",{staticClass:"language-bash extra-class"},[s("pre",{pre:!0,attrs:{class:"language-bash"}},[s("code",[a._v("conda "),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("install")]),a._v(" -c conda-forge -c bioconda magpurify2\n")])])])])],1),a._v(" "),s("p",[a._v("If you choose install MAGpurify2 via "),s("code",[a._v("pip")]),a._v(", make sure that you have also installed the third-party dependencies: "),s("a",{attrs:{href:"https://github.com/hyattpd/Prodigal",target:"_blank",rel:"noopener noreferrer"}},[a._v("Prodigal"),s("OutboundLink")],1),a._v(" and "),s("a",{attrs:{href:"https://github.com/soedinglab/MMseqs2",target:"_blank",rel:"noopener noreferrer"}},[a._v("MMseqs2"),s("OutboundLink")],1),a._v(". The "),s("code",[a._v("conda")]),a._v(" install method will automatically download and install these software.")]),a._v(" "),s("h2",{attrs:{id:"download-test-data-and-database"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#download-test-data-and-database"}},[a._v("#")]),a._v(" Download test data and database")]),a._v(" "),s("ul",[s("li",[a._v("Download test data:")])]),a._v(" "),s("div",{staticClass:"language-bash extra-class"},[s("pre",{pre:!0,attrs:{class:"language-bash"}},[s("code",[s("span",{pre:!0,attrs:{class:"token assign-left variable"}},[a._v("fileId")]),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"1-Gf-FsVIcARqrUb-LHS_FZlb-sCGAmfo"')]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token assign-left variable"}},[a._v("fileName")]),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"magpurify2_test_data.tar.gz"')]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("curl")]),a._v(" -sc /tmp/cookie "),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"https://drive.google.com/uc?export=download&id='),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${fileId}")]),a._v('"')]),a._v(" "),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v(">")]),a._v(" /dev/null\n"),s("span",{pre:!0,attrs:{class:"token assign-left variable"}},[a._v("code")]),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"'),s("span",{pre:!0,attrs:{class:"token variable"}},[s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("$(")]),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("awk")]),a._v(" "),s("span",{pre:!0,attrs:{class:"token string"}},[a._v("'/_warning_/ {print "),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("$NF")]),a._v("}'")]),a._v(" /tmp/cookie"),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v(")")])]),a._v('"')]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("curl")]),a._v(" -Lb /tmp/cookie "),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"https://drive.google.com/uc?export=download&confirm='),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${code}")]),a._v("&id="),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${fileId}")]),a._v('"')]),a._v(" -o "),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${fileName}")]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("tar")]),a._v(" zxfv magpurify2_test_data.tar.gz\n")])])]),s("ul",[s("li",[a._v("Download database:")])]),a._v(" "),s("div",{staticClass:"language-bash extra-class"},[s("pre",{pre:!0,attrs:{class:"language-bash"}},[s("code",[s("span",{pre:!0,attrs:{class:"token assign-left variable"}},[a._v("fileId")]),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"1ooWiR3LplBy5GsY5wZ7o6dwswiCWVvmi"')]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token assign-left variable"}},[a._v("fileName")]),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"magpurify2DB.v1.0.tar.gz"')]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("curl")]),a._v(" -sc /tmp/cookie "),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"https://drive.google.com/uc?export=download&id='),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${fileId}")]),a._v('"')]),a._v(" "),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v(">")]),a._v(" /dev/null\n"),s("span",{pre:!0,attrs:{class:"token assign-left variable"}},[a._v("code")]),s("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"'),s("span",{pre:!0,attrs:{class:"token variable"}},[s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("$(")]),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("awk")]),a._v(" "),s("span",{pre:!0,attrs:{class:"token string"}},[a._v("'/_warning_/ {print "),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("$NF")]),a._v("}'")]),a._v(" /tmp/cookie"),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v(")")])]),a._v('"')]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("curl")]),a._v(" -Lb /tmp/cookie "),s("span",{pre:!0,attrs:{class:"token string"}},[a._v('"https://drive.google.com/uc?export=download&confirm='),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${code}")]),a._v("&id="),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${fileId}")]),a._v('"')]),a._v(" -o "),s("span",{pre:!0,attrs:{class:"token variable"}},[a._v("${fileName}")]),a._v("\n"),s("span",{pre:!0,attrs:{class:"token function"}},[a._v("tar")]),a._v(" zxfv magpurify2DB.v1.0.tar.gz\n")])])]),s("h2",{attrs:{id:"docker"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#docker"}},[a._v("#")]),a._v(" Docker")]),a._v(" "),s("p",[a._v("Foo.")])],1)}),[],!1,null,null,null);t.default=r.exports}}]);