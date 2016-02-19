import qbs 1.0

Product {
    type: "application"
    consoleApplication: true
    name : "manual-tests"
    files : [
        "henyeygreensteinmanual_tests.cpp",
        "henyeygreensteinvector_tests.cpp",
        "main.cpp",
        "schema/NeuroML_v2beta3.cpp",
        "schema/NeuroML_v2beta3.h",
        "voxelizer.cpp",
    ]
    cpp.cxxLanguageVersion: "c++14"
    cpp.cppFlags: ["-fopenmp"]
    cpp.linkerFlags: ["-fopenmp"]
    cpp.dynamicLibraries: ["xerces-c", "hdf5"]
    cpp.libraryPaths: ["/usr/lib/x86_64-linux-gnu/hdf5/serial"]
    cpp.includePaths: ["/usr/lib/x86_64-linux-gnu/hdf5/serial/include"]

    Depends { name: "Qt.core" }
    Depends { name: "Qt.xml" }

    Depends { name: "cpp" }
    Depends { name: "photonflow-lib" }
    Depends { name: "conanbuildinfo" }
}
