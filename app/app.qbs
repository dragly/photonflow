import qbs 1.0

Product {
    type: "application"
    name : "photonflow-app"
    files : [
        "main.cpp",
        "qml.qrc",
        "visual/renderview.cpp",
        "visual/renderview.h",
    ]
    cpp.cxxFlags: ["-fopenmp"]
    cpp.linkerFlags: ["-fopenmp"]
    cpp.cxxLanguageVersion: "c++14"
    Depends {
        name: "photonflow-lib"
    }
    Depends {
        name: "Qt"
        submodules: [
            "core",
            "gui",
            "qml",
            "quick",
            "widgets"
        ]
    }
}
