import SimVis 1.0
import SimVis.ShaderNodes 1.0
import Qt3D.Core 2.0
import Qt3D.Logic 2.0
import Qt3D.Render 2.0
import QtQuick.Scene3D 2.0
import QtQuick 2.5 as QQ2
import Photonflow 1.0

Scene3D {
    id: root

    property alias camera: visualizer.camera
    property var currentEntity
    property var neurons: []

    aspects: ["input", "render", "logic"]
    focus: true

    function addNeuron(path) {
        var neuronComponent = Qt.createComponent("Neuron.qml")
//        var properties = {
//            transform: {
//                translation: Qt.vector3d(0, 0, 0)
//            }
//        }
        if(neuronComponent.status !== QQ2.Component.Ready) {
            console.log("Could not create neuron")
            throw(neuronComponent.errorString())
        }

        var neuron = neuronComponent.createObject(visualizer)
        neuron.parent = visualizer
        neuron.pressed.connect(function() {
            console.log("Setting current entity to", neuron)
            for(var i in neurons) {
                var otherNeuron = neurons[i]
                if(otherNeuron !== neuron) {
                    otherNeuron.selected = false
                }
            }
            neuron.selected = true
            currentEntity = neuron
        })
        neurons.push(neuron)
        neuron.pressed()
    }

    Visualizer {
        id: visualizer

        camera.aspectRatio: root.width / root.height
        camera.position: Qt.vector3d(8, 8, 8)

        CameraController {
            camera: visualizer.camera
        }

        EntityController {
            id: entityController
            camera: visualizer.camera
            entity: root.currentEntity
        }
    }
}
