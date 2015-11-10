import qbs 1.0

Product {
    type: "application"
    consoleApplication: true
    name : "photonflow-tests"
    files : [ "main.cpp" ]
    cpp.cxxLanguageVersion: "c++14"
    Depends { name: "cpp" }
    Depends { name: "photonflow-lib" }
}
