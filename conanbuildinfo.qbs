import qbs 1.0

Project {
    Product {
        name: "conanbuildinfo"
        Export {
            Depends { name: "cpp" }
            cpp.includePaths: ["/home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed/include"]
            cpp.libraryPaths: ["/home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed/lib"]
            cpp.systemIncludePaths: ["/home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed/bin"]
            cpp.dynamicLibraries: []
            cpp.defines: []
            cpp.cppFlags: []
            cpp.cFlags: []
            cpp.linkerFlags: []
        }
    }

    Product {
        name: "Catch"
        Export {
            Depends { name: "cpp" }
            cpp.includePaths: ["/home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed/include"]
            cpp.libraryPaths: ["/home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed/lib"]
            cpp.systemIncludePaths: ["/home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed/bin"]
            cpp.dynamicLibraries: []
            cpp.defines: []
            cpp.cppFlags: []
            cpp.cFlags: []
            cpp.linkerFlags: []
        }
    }
    // Catch root path: /home/svenni/.conan/data/Catch/1.3.2/dragly/master/package/0692fb2bd888ba708ca65670557c56d2e16851ed
}
