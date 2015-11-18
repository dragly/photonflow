import qbs 1.0

Product {
    type: "application"
    consoleApplication: true
    name : "auto-tests"
    files : [
        "henyeygreenstein_tests.cpp",
        "henyeygreensteindistribution_tests.cpp",
        "main.cpp",
    ]
    cpp.cxxLanguageVersion: "c++14"
    Depends { name: "cpp" }
    Depends { name: "photonflow-lib" }
}
