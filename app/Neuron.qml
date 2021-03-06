import QtQuick.Scene3D 2.0
import QtQuick 2.0 as QQ2

import Qt3D.Core 2.0
import Qt3D.Logic 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.0

import SimVis 1.0
import ShaderNodes 1.0
import Photonflow 1.0

Entity {
    id: root
    signal pressed(var pick)

    property bool selected: false
    property alias transform: transform_
    property alias simulator: simulator_

    property var lights: [light1, light2, light3, light4]

    Light {
        id: light1
        position: Qt.vector3d(-100, 100, -100)
        strength: 0.4
        attenuation: 0.0
    }
    Light {
        id: light2
        position: Qt.vector3d(-100, 100, 100)
        strength: 0.4
        attenuation: 0.0
    }
    Light {
        id: light3
        position: Qt.vector3d(100, 100, 100)
        strength: 0.4
        attenuation: 0.0
    }
    Light {
        id: light4
        position: Qt.vector3d(100, 100, -100)
        strength: 0.4
        attenuation: 0.0
    }

    components: [
        Transform {
            id: transform_
        },
        SphereMesh {
            radius: 0.1
        },
        ShaderBuilderMaterial {
            fragmentColor: StandardMaterial {
                color: root.selected ? "#FF0000" : "#0000FF"
                lights: root.lights //TODO fix lights with new api
            }
        },
        ObjectPicker {
            hoverEnabled: true
            onPressed: {
                root.pressed(pick)
            }
        }
    ]

    NeuronSimulator {
        id: simulator_
    }
    Cylinders {
        cylinderData: simulator_.cylinderData
        fragmentColor: StandardMaterial {
            color: root.selected ? "#F7C7C7" : "#E7E7E7"
            lights: root.lights
        }
    }
}
