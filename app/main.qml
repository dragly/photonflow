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
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(0, 0, 0)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(1, 0, 0)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(0, 1, 0)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(1, 1, 0)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(0, 0, 1)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(1, 0, 1)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(0, 1, 1)})
        builderScene.addNeuron("...", {"transform.translation": Qt.vector3d(1, 1, 1)})
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

        camera: builderScene.simulatorCamera
    }

    RowLayout {
        anchors.fill: parent

        Rectangle {
            Layout.fillHeight: true
            Layout.fillWidth: true

            width: 1

            ColumnLayout {
                anchors.fill: parent

                Text {
                    text: "Render time: " + simulator.renderTime
                }

                Text {
                    text: "Mode: " + builderScene.mode
                }

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
                    maximumValue: 1000.0
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
                    stepSize: 0.001
                    minimumValue: 0.94
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
        }

        ImageViewer {
            Layout.fillHeight: true
            Layout.fillWidth: true
            width: 3
            image: simulator.image
        }
    }
}

