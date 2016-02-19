import qbs 1.0

Project {
    references: [
        "app/app.qbs",
        "lib/lib.qbs",
        "tests/tests.qbs",
        "conanbuildinfo.qbs"
    ]
}
