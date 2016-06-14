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
    property Entity currentEntity
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
        neuron.pressed.connect(function(pick) {
            selectEntity(neuron)
            mouseHandler.pickOffset = pick.worldIntersection.minus(neuron.transform.translation)
            mouseHandler.pickPosition = pick.worldIntersection
            mouseHandler.pickRegistered = true
            console.log("Neuron pressed")
        })
        entities.push(neuron)
        neurons.push(neuron)
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
            id: mouseHandler

            property vector3d pickOffset
            property vector3d pickPosition
            property vector3d previousPosition
            property bool pickRegistered: false
            property bool pressRegistered: false

            sourceDevice: mouseDevice

            function intersection(mouse) {
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
                var x0 = pickPosition
                var d = -x0.dotProduct(n)
                var p0 = camera.position
                var t = -(p0.dotProduct(n) + d) / (v.dotProduct(n))

                var p = p0.plus(v.times(t))

                return p
            }

            onPressed: {
                console.log("Pressed", pickRegistered)
                if(pickRegistered) {
                    previousPosition = intersection(mouse)
                    pressRegistered = true
                }
            }

            onPositionChanged: {
                if(!currentEntity || !pickRegistered) {
                    return
                }

                var position = intersection(mouse)
                var diff = position.minus(previousPosition)

                if(pressRegistered) {
                    currentEntity.transform.translation = currentEntity.transform.translation.plus(diff)
                }

                previousPosition = position
                pressRegistered = true
            }

            onReleased: {
                pickRegistered = false
                pressRegistered = false
            }
        }

        CameraController {
            enabled: !mouseHandler.pickRegistered
            mouseSourceDevice: mouseDevice
            keyboardSourceDevice: keyboardDevice
            camera: visualizer.camera
        }

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
