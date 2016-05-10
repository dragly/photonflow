import qbs 1.0

Project {
    references: [
        "app/app.qbs",
        "lib/lib.qbs",
        "tests/tests.qbs",
        "tools/tools.qbs",

        "conanbuildinfo.qbs",
        "vendor/SimVis/package.qbs"
    ]
}

