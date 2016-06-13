import SimVis 1.0
import SimVis.ShaderNodes 1.0
import Qt3D.Core 2.0
import Qt3D.Logic 2.0
import Qt3D.Input 2.0
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
        simulatorCamera_
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

    function addNeuron(path, properties) {
        var neuronComponent = Qt.createComponent("Neuron.qml")
        if(neuronComponent.status !== QQ2.Component.Ready) {
            console.log("Could not create neuron")
            throw(neuronComponent.errorString())
        }

        if(!properties) {
            properties = {}
        }

        var neuron = neuronComponent.createObject(visualizer, properties)
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
        camera.position: Qt.vector3d(0, 0, 5)

        KeyboardDevice {
            id: keyboardDevice
        }

        MouseDevice {
            id: mouseDevice
            sensitivity: 0.01
        }

        MouseHandler {
            sourceDevice: mouseDevice

            onPositionChanged: {
                var camera = visualizer.camera
                var viewProjection = camera.projectionMatrix.times(camera.viewMatrix)
                var viewProjectionInverse = viewProjection.inverted()
                var screenX = 2.0 * mouse.x / root.width - 1.0
                var screenY = 1.0 - 2.0 * mouse.y / root.height
                var screenVector = Qt.vector4d(screenX, screenY, -1.0, 1.0)
                var worldVectorTemp = viewProjectionInverse.times(screenVector)
                var worldVector = worldVectorTemp.times(1.0 / worldVectorTemp.w)
                var ray = worldVector.toVector3d().minus(camera.position)

                // ray-plane intersection
                var n = camera.viewVector
                var v = ray
                var x0 = currentEntity.transform.translation
                var d = -x0.dotProduct(n)
                var p0 = camera.position
                var t = -(p0.dotProduct(n) + d) / (v.dotProduct(n))

                var p = p0.plus(v.times(t))

                console.log("p", p)

                currentEntity.transform.translation = p

                // TODO use difference since last position instead
            }
        }

//        CameraController {
//            mouseSourceDevice: mouseDevice
//            keyboardSourceDevice: keyboardDevice
//            camera: visualizer.camera
//        }

        EntityController {
            id: entityController
            mouseSourceDevice: mouseDevice
            keyboardSourceDevice: keyboardDevice
            camera: visualizer.camera
            entity: root.currentEntity
        }

        Entity {
            id: simulatorCamera_
            property bool selected: false
            property alias transform: cameraBoxTransform
            components: [
                Mesh {
                    source: "meshes/cameracone.obj"
                },
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

        Test {
        }
    }
}
