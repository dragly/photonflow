import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import QtQuick.Layouts 1.1
import Qt.labs.settings 1.0

import Photonflow 1.0

ApplicationWindow {
    visible: true
    width: 1280
    height: 1024
    title: qsTr("Photonflow")

    Settings {
        property alias emissionFactor: simulator.emissionFactor
        property alias absorptionCoefficient: simulator.absorptionCoefficient
        property alias henyeyGreensteinFactor: simulator.henyeyGreensteinFactor
        property alias scatteringCoefficient: simulator.scatteringCoefficient
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

                Text {
                    text: "Emission: " + simulator.emissionFactor.toFixed(2)
                }

                Slider {
                    id: emissionSlider
                    Layout.fillWidth: true
                    value: simulator.emissionFactor
                    minimumValue: 0.01
                    maximumValue: 1.0

                    Binding {
                        target: simulator
                        property: "emissionFactor"
                        value: emissionSlider.value
                    }
                }

                Text {
                    text: "Scattering: " + simulator.scatteringCoefficient.toFixed(2)
                }

                Slider {
                    id: scatteringSlider
                    Layout.fillWidth: true
                    value: simulator.scatteringCoefficient
                    minimumValue: 0.0
                    maximumValue: 1.0

                    Binding {
                        target: simulator
                        property: "scatteringCoefficient"
                        value: scatteringSlider.value
                    }
                }

                Text {
                    text: "Absorption: " + simulator.absorptionCoefficient.toFixed(2)
                }

                Slider {
                    id: absorptionSlider
                    Layout.fillWidth: true
                    value: simulator.absorptionCoefficient
                    minimumValue: 0.0
                    maximumValue: 1.0

                    Binding {
                        target: simulator
                        property: "absorptionCoefficient"
                        value: absorptionSlider.value
                    }
                }

                Text {
                    text: "Henyey greenstein factor: " + simulator.henyeyGreensteinFactor.toFixed(2)
                }

                Slider {
                    id: henyeyGreensteinFactorSlider
                    Layout.fillWidth: true
                    value: simulator.henyeyGreensteinFactor
                    minimumValue: 0.0
                    maximumValue: 1.0

                    Binding {
                        target: simulator
                        property: "henyeyGreensteinFactor"
                        value: henyeyGreensteinFactorSlider.value
                    }
                }

                Item {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                }
            }
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

