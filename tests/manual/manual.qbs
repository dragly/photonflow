import qbs 1.0

Product {
    type: "application"
    consoleApplication: true
    name : "manual-tests"
    files : [
        "henyeygreensteinmanual_tests.cpp",
        "henyeygreensteinvector_tests.cpp",
        "main.cpp",
    ]
    cpp.cxxLanguageVersion: "c++14"
    Depends { name: "cpp" }
    Depends { name: "photonflow-lib" }
}
