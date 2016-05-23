import SimVis 1.0
import SimVis.ShaderNodes 1.0
import Qt3D.Core 2.0
import Qt3D.Logic 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.0
import QtQuick.Scene3D 2.0
import QtQuick 2.5 as QQ2
import Photonflow 1.0

Scene3D {
    id: root

    property alias simulatorCamera: simulatorCamera_
    property alias mode: entityController.mode
    property alias viewportCamera: visualizer.camera
    property var currentEntity
    property var neurons: []
    property var entities: [
        simulatorCamera_,
        boundingBox
    ]

    aspects: ["input", "render", "logic"]
    focus: true

    function selectEntity(entity) {
        console.log("Setting current entity to", entity)
        for(var i in entities) {
            var otherEntity = entities[i]
            if(otherEntity !== entity) {
                otherEntity.selected = false
            }
        }
        entity.selected = true
        currentEntity = entity
    }

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
            selectEntity(neuron)
        })
        entities.push(neuron)
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

        Entity {
            id: boundingBox
            property alias transform: boundingBoxTransform
            property bool selected: false
            components: [
                CuboidMesh {},
                PhongMaterial {
                    diffuse: boundingBox.selected ? "red" : "lightblue"
                },
                Transform {
                    id: boundingBoxTransform
                    scale3D: Qt.vector3d(1, 2, 2)
                },
                ObjectPicker {
                    hoverEnabled: true
                    onPressed: root.selectEntity(boundingBox)
                }
            ]
        }

        Entity {
            id: simulatorCamera_
            property bool selected: false
            property alias transform: cameraBoxTransform
            components: [
                CuboidMesh {},
                PhongMaterial {
                    diffuse: simulatorCamera_.selected ? "red" : "lightblue"
                },
                Transform {
                    id: cameraBoxTransform
                    translation: Qt.vector3d(0, 0, -10)
                },
                ObjectPicker {
                    hoverEnabled: true
                    onPressed: root.selectEntity(simulatorCamera_)
                }
            ]
        }
    }
}
