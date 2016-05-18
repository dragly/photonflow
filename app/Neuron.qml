import QtQuick.Scene3D 2.0

import Qt3D.Core 2.0
import Qt3D.Logic 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.0

import SimVis 1.0
import SimVis.ShaderNodes 1.0
import Photonflow 1.0

Entity {
    id: root
    signal pressed

    property bool selected: false
    property alias transform: transform_
    property alias simulator: simulator_

    components: [
        Transform {
            id: transform_
        },
        SphereMesh {
            radius: 0.1
        },
        PhongMaterial {
            diffuse: root.selected ? "red" : "lightblue"
        },
        ObjectPicker {
            hoverEnabled: true
            onPressed: root.pressed()
        }
    ]
    NeuronSimulator {
        id: simulator_
    }
    Cylinders {
        cylinderData: simulator_.cylinderData
    }
}
