import qbs 1.0

Product {
    id: libRoot
    property string includesInstallPath: "include/photonflow"

    type: "staticlibrary"
    name: "photonflow-lib"

    files: [
        "armadillo_includer.h",
        "cameras/perspective.cpp",
        "cameras/perspective.h",
        "core/camera.h",
        "core/camera.cpp",
        "core/common.h",
        "core/diffgeom.cpp",
        "core/diffgeom.h",
        "core/error.cpp",
        "core/error.h",
        "core/film.cpp",
        "core/film.h",
        "core/filter.cpp",
        "core/filter.h",
        "core/geometry.cpp",
        "core/geometry.h",
        "core/heyneygreenstein.h",
        "core/integrator.cpp",
        "core/integrator.h",
        "core/intersection.cpp",
        "core/intersection.h",
        "core/memory.cpp",
        "core/memory.h",
        "core/montecarlo.cpp",
        "core/montecarlo.h",
        "core/quaternion.cpp",
        "core/quaternion.h",
        "core/randomnumbergenerator.cpp",
        "core/randomnumbergenerator.h",
        "core/renderer.cpp",
        "core/renderer.h",
        "core/sampler.cpp",
        "core/sampler.h",
        "core/scene.cpp",
        "core/scene.h",
        "core/shape.cpp",
        "core/shape.h",
        "core/spectrum.cpp",
        "core/spectrum.h",
        "core/transform.cpp",
        "core/transform.h",
        "core/units.cpp",
        "core/units.h",
        "core/volume.cpp",
        "core/volume.h",
        "film/image.cpp",
        "film/image.h",
        "filters/box.cpp",
        "filters/box.h",
        "filters/gaussian.cpp",
        "filters/gaussian.h",
        "filters/mitchell.cpp",
        "filters/mitchell.h",
        "filters/sinc.cpp",
        "filters/sinc.h",
        "filters/triangle.cpp",
        "filters/triangle.h",
        "geometry/bbox.cpp",
        "geometry/bbox.h",
        "geometry/cylinderfrustum.cpp",
        "geometry/cylinderfrustum.h",
        "geometry/normal.cpp",
        "geometry/normal.h",
        "geometry/point3d.cpp",
        "geometry/point3d.h",
        "geometry/ray.cpp",
        "geometry/ray.h",
        "geometry/raydifferential.cpp",
        "geometry/raydifferential.h",
        "geometry/rectangle.cpp",
        "geometry/rectangle.h",
        "geometry/vector3d.cpp",
        "geometry/vector3d.h",
        "io/neuromlreader.cpp",
        "io/neuromlreader.h",
        "io/voxelizer.cpp",
        "io/voxelizer.h",
        "samplers/random.cpp",
        "samplers/random.h",
        "stdafx.h",
        "tools/vsds.py",
        "volumes/volumegrid.cpp",
        "volumes/volumegrid.h",
    ]
    cpp.includePaths: ["."]
    cpp.cxxLanguageVersion: "c++14"
    Depends { name: 'cpp' }
    Depends {
        name: "Qt"
        submodules: ["core", "gui", "xml"]
    }

    Export {
        Depends { name: "cpp" }
        cpp.includePaths: [
            "."
        ]
    }
}
