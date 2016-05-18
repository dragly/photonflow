import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import QtQuick.Layouts 1.1
import Qt.labs.settings 1.0

import Photonflow 1.0

ApplicationWindow {
    id: root

    visible: true
    width: 1920
    height: 1024
    title: qsTr("Photonflow")

    Component.onCompleted: {
        builderScene.addNeuron("...")
        builderScene.addNeuron("...")
    }

    function voxelize() {
        var neuronSimulators = []
        for(var i in builderScene.neurons) {
            var neuron = builderScene.neurons[i]
            neuronSimulators.push({simulator: neuron.simulator, transform: neuron.transform})
        }
        simulator.voxelize(neuronSimulators)
    }

    Settings {
        property alias emissionFactor: simulator.emissionFactor
        property alias absorptionCoefficient: simulator.absorptionCoefficient
        property alias henyeyGreensteinFactor: simulator.henyeyGreensteinFactor
        property alias scatteringCoefficient: simulator.scatteringCoefficient
        property alias windowX: root.x
        property alias windowY: root.y
        property alias windowWidth: root.width
        property alias windowHeight: root.height
    }

    PhotonflowSimulator {
        id: simulator
    }

    RowLayout {
        anchors.fill: parent

        Rectangle {
            Layout.fillHeight: true
            Layout.fillWidth: true

            width: 1

            ColumnLayout {
                anchors.fill: parent

                Button {
                    text: simulator.running ? "Stop" : "Start"
                    onClicked: {
                        simulator.running = !simulator.running
                    }
                }

                Button {
                    text: "Clear"
                    onClicked: {
                        simulator.clear()
                    }
                }

                Button {
                    text: "Create neuron"
                    onClicked: {
                        builderScene.addNeuron("...")
                    }
                }

                Button {
                    text: "Voxelize"
                    onClicked: {
                        voxelize()
                    }
                }

                BoundSlider {
                    id: emissionSlider
                    Layout.fillWidth: true
                    text: "Emission"
                    minimumValue: 0.01
                    maximumValue: 1.0
                    target: simulator
                    property: "emissionFactor"
                }

                BoundSlider {
                    id: scatteringSlider
                    Layout.fillWidth: true
                    text: "Scattering"
                    minimumValue: 0.0
                    maximumValue: 1.0
                    target: simulator
                    property: "scatteringCoefficient"
                }

                BoundSlider {
                    id: absorptionSlider
                    text: "Absorption"
                    Layout.fillWidth: true
                    minimumValue: 0.0
                    maximumValue: 1.0
                    target: simulator
                    property: "absorptionCoefficient"
                }

                BoundSlider {
                    id: henyeyGreensteinFactorSlider
                    target: simulator
                    property: "henyeyGreensteinFactor"
                    text: "Henyey Greenstein factor"
                    minimumValue: 0.0
                    maximumValue: 1.0
                    Layout.fillWidth: true
                }

                Item {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                }
            }
        }

        Item {
            Layout.fillWidth: true
            Layout.fillHeight: true
            width: 3
            BuilderScene {
                id: builderScene
                anchors.fill: parent
            }

//            MouseArea {
//                anchors.fill: parent
//                propagateComposedEvents: true
//                onPressed: {
//                    console.log("Pressed")
////                    mouse.accepted = false
//                }
//                onReleased: {
//                    console.log("Released")
//                    mouse.accepted = false
//                }
//            }

//            MouseArea {
//                property point previousPosition
//                property real dragSpeed: 0.04
//                property alias camera: builderScene.camera
//                property alias entity: builderScene.currentEntity

//                propagateComposedEvents: true
//                acceptedButtons: Qt.LeftButton | Qt.RightButton
//                anchors.fill: parent
//                onPressed: {
//                    previousPosition = Qt.point(mouse.x, mouse.y)
//                }
//                onPositionChanged: {
//                    var currentPosition = Qt.point(mouse.x, mouse.y)
//                    var diff = Qt.vector2d(previousPosition.x - currentPosition.x, previousPosition.y - currentPosition.y)

//                    diff = diff.times(0.5)

//                    if(mouse.buttons & Qt.LeftButton) {
//                        camera.panAboutViewCenter(diff.x, camera.upVector)
//                        camera.tiltAboutViewCenter(-diff.y)
//                    } else if(mouse.buttons & Qt.RightButton) {
//                        var rightVector = camera.viewVector.crossProduct(camera.upVector).normalized()
//                        var upVector = camera.upVector.normalized()
//                        var direction = rightVector.times(-diff.x).plus(upVector.times(diff.y))
//                        entity.transform.translation = entity.transform.translation.plus(direction.times(dragSpeed))
//                    }
//                    previousPosition = Qt.point(mouse.x, mouse.y)
//                }
//                onReleased: {
//                    mouse.accepted = false
//                }
//                onClicked: {
//                    mouse.accepted = false
//                }
//                onWheel: {
//                    camera.fieldOfView = Math.max(10.0, Math.min(160.0, camera.fieldOfView - wheel.angleDelta.y * 0.1))
//                }
//            }
        }

        ImageViewer {
            Layout.fillHeight: true
            Layout.fillWidth: true
            width: 3
            image: simulator.image
        }

    }

//    Timer {
//        running: true
//        repeat: true
//        interval: 100
//        onTriggered: {
//            simulator.step()
//        }
//    }

//    Visualizer {
//        anchors.fill: parent
//        simulator: TestSimulator {
//            id: simulator
//        }
//    }

//    Text {
//        color: "white"
//        text: Math.round(timer.totalTime, 1)
//    }

//    MouseArea {
//        anchors.fill: parent
//        onClicked: timer.running = !timer.running
//    }

//    Timer {
//        id: timer
//        property real lastTime: Date.now()
//        property real totalTime: 0
//        running: true
//        repeat: true
//        interval: 16

//        onTriggered: {
//            console.log("Triggered!")
//            totalTime += Date.now() - lastTime
//            lastTime = Date.now()
//            renderView.requestIntegrate()
//        }
//    }
}

