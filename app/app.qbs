import qbs

Product {
    property string hdf5LibraryName: "hdf5_serial"
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
    cpp.dynamicLibraries: [hdf5LibraryName]
    cpp.cxxLanguageVersion: "c++14"
    Depends {
        name: "photonflow-lib"
    }
    Depends {
        name: "Qt"
        submodules: [
            "core",
            "gui",
            "opengl",
            "qml",
            "quick",
            "widgets",
            "3dcore",
            "3drender",
            "3dinput",
            "3dquick"
        ]
    }
}

