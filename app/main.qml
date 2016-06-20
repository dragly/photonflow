import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import QtQuick.Layouts 1.1
import QtGraphicalEffects 1.0
import Qt.labs.settings 1.0

import Photonflow 1.0

ApplicationWindow {
    id: root

    visible: true
    width: 1920
    height: 1024
    title: qsTr("Photonflow")

    Component.onCompleted: {
        var dummyTransform = Qt.createQmlObject("import Qt3D.Core 2.0; Transform {}", root)

        for(var i = 0; i < 2; i++) {
            for(var j = 0; j < 2; j++) {
                for(var k = 0; k < 2; k++) {
                    builderScene.addNeuron("...", {
                                               "transform.translation": Qt.vector3d(i - 1 + Math.random(), j - 2 + Math.random(), k - 1 + Math.random()),
                                               "transform.rotation": dummyTransform.fromAxisAndAngle(Qt.vector3d(0, 1, 0), Math.random() * 360)
                                           })
                }
            }
        }

        //        builderScene.addNeuron("...", {})

        dummyTransform.destroy()
    }

    function rebuild() {
        simulator.running = true
        modeMenu.currentIndex = 1
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
        running: false
    }

    RowLayout {
        anchors.fill: parent
        spacing: 0

        Rectangle {
            Layout.preferredWidth: 96
            Layout.fillHeight: true

            LinearGradient {
                anchors.fill: parent
                start: Qt.point(0, 0)
                end: Qt.point(width, 10)
                gradient: Gradient {
                    GradientStop {
                        color: "#222"
                        position: 0.0
                    }
                    GradientStop {
                        color: "#1a1a1a"
                        position: 1.0
                    }
                }
            }

            Image {
                id: photonflowLogo
                anchors {
                    top: parent.top
                    left: parent.left
                    right: parent.right
                    margins: 12
                }

                height: width
                source: "qrc:/images/logo/photonflow-logo_48dp.png"
                fillMode: Image.PreserveAspectFit
            }

            Text {
                id: logoText
                anchors {
                    horizontalCenter: parent.horizontalCenter
                    top: photonflowLogo.bottom
                    topMargin: 8
                }

                text: "photonflow\n0.1"
                color: "#888"
                font.pixelSize: 12
                horizontalAlignment: Text.AlignHCenter
            }

            Column {
                id: modeMenu

                property int currentIndex: 0

                anchors {
                    verticalCenter: parent.verticalCenter
                    left: parent.left
                    right: parent.right
                    margins: 16
                }
                spacing: 24

                Repeater {
                    model: ListModel {
                        ListElement {
                            image: "qrc:/images/ic_brush_white_48dp.png"
                            name: "Model"
                        }
                        ListElement {
                            image: "qrc:/images/ic_image_white_48dp.png"
                            name: "Render"
                        }
                        ListElement {
                            image: "qrc:/images/ic_timeline_white_48dp.png"
                            name: "Analyze"
                        }
                        ListElement {
                            image: "qrc:/images/ic_help_white_48dp.png"
                            name: "Help"
                        }
                    }

                    ToggleButton {
                        anchors {
                            left: parent.left
                            right: parent.right
                        }
                        toggled: modeMenu.currentIndex === index
                        source: image
                        text: name
                        onClicked: {
                            modeMenu.currentIndex = index
                        }
                    }
                }
            }

            Column {
                anchors {
                    bottom: parent.bottom
                    left: parent.left
                    right: parent.right
                    margins: 16
                }
                spacing: 24

                ToggleButton {
                    anchors {
                        left: parent.left
                        right: parent.right
                        margins: 8
                    }
                    source: "qrc:/images/ic_play_arrow_white_48dp.png"
                    toggled: simulator.running
                    text: "Simulate"
                    onClicked: {
                        simulator.running = !simulator.running
                    }
                }
                ToggleButton {
                    anchors {
                        left: parent.left
                        right: parent.right
                        margins: 8
                    }
                    source: "qrc:/images/ic_refresh_white_48dp.png"
                    text: "Rebuild"
                    onClicked: {
                        rebuild()
                    }
                }
            }
        }

        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true

            color: "#333"

            Text {
                id: titleText
                anchors {
                    horizontalCenter: viewport.horizontalCenter
                    top: parent.top
                    margins: 16
                }

                text: "Project 1"
                color: "#cccccc"
                font.pixelSize: 32
                font.weight: Font.Light
            }

            Item {
                id: viewport
                anchors {
                    top: titleText.bottom
                    bottom: parent.bottom
                    left: parent.left
                    right: toolbar.left
                    margins: 16
                }

                BuilderScene {
                    id: builderScene
                    anchors.fill: parent
                    visible: modeMenu.currentIndex === 0
                    clearColor: "#444"
                }

                ImageViewer {
                    visible: modeMenu.currentIndex === 1
                    anchors.fill: parent
                    image: simulator.image
                }
            }

            Item {
                id: toolbar

                anchors {
                    top: titleText.bottom
                    right: parent.right
                    bottom: parent.bottom
                    margins: 16
                }

                width: 360

                Row {
                    id: tabMenu

                    property int currentIndex: 0
                    property string currentIdentifier: tabRepeater.model.get(currentIndex).identifier

                    anchors {
                        left: parent.left
                        right: parent.right
                    }
                    height: 48
                    spacing: 8

                    Repeater {
                        id: tabRepeater
                        anchors.fill: parent

                        model: ListModel {
                            ListElement {
                                image: "qrc:/images/ic_device_hub_white_48dp.png"
                                identifier: "neuron"
                            }
                            ListElement {
                                image: "qrc:/images/ic_blur_on_white_48dp.png"
                                identifier: "volume"
                            }
                            ListElement {
                                image: "qrc:/images/ic_camera_white_48dp.png"
                                identifier: "microscope"
                            }
                            ListElement {
                                image: "qrc:/images/ic_wb_sunny_white_48dp.png"
                                identifier: "fluorescence"
                            }
                            ListElement {
                                image: "qrc:/images/ic_image_white_48dp.png"
                                identifier: "render"
                            }
                        }

                        ToggleImage {
                            width: height
                            height: tabMenu.height - tabMenu.spacing
                            toggled: tabMenu.currentIndex === index
                            source: image
                            MouseArea {
                                anchors.fill: parent
                                onClicked: tabMenu.currentIndex = index
                            }
                        }
                    }
                }

                Rectangle {
                    anchors {
                        top: tabMenu.bottom
                        topMargin: 8
                        left: parent.left
                        right: parent.right
                        bottom: parent.bottom
                    }

                    color: "#222222"
                    radius: 8

                    Item {
                        anchors {
                            fill: parent
                            margins: 8
                        }

                        Column {
                            anchors.fill: parent
                            visible: tabMenu.currentIdentifier === "render"
                            spacing: 8

                            Text {
                                color: "#ccc"
                                text: "Resolution:"
                            }

                            RowLayout {
                                anchors {
                                    left: parent.left
                                    right: parent.right
                                }
                                height: widthTextInput.height

                                TextField {
                                    id: widthTextInput
                                    Layout.fillWidth: true
                                    text: simulator.resolution.width

                                    Binding {
                                        target: simulator
                                        property: "resolution.width"
                                        value: widthTextInput.text
                                    }
                                }

                                Text {
                                    text: "x"
                                }

                                TextField {
                                    id: heightTextInput
                                    Layout.fillWidth: true
                                    text: simulator.resolution.height

                                    Binding {
                                        target: simulator
                                        property: "resolution.height"
                                        value: heightTextInput.text
                                    }
                                }
                            }

                            Text {
                                text: "Preview:"
                                color: "#ccc"
                            }

                            ImageViewer {
                                anchors {
                                    left: parent.left
                                    right: parent.right
                                }

                                height: width * 9 / 16

                                image: simulator.image
                            }

                            Text {
                                text: "Samples: " + (simulator.completedSampleCount / 1e6).toFixed(1) + "M"
                                color: "#ccc"
                            }

                            Text {
                                text: "Render time: " + simulator.renderTime
                                color: "#ccc"
                            }
                        }

                        Column {
                            anchors.fill: parent
                            visible: tabMenu.currentIdentifier === "volume"

                            BoundSlider {
                                id: scatteringSlider

                                text: "Scattering"
                                minimumValue: 0.0
                                maximumValue: 1.0
                                target: simulator
                                property: "scatteringCoefficient"
                            }

                            BoundSlider {
                                id: absorptionSlider
                                text: "Absorption"
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
                            }
                        }

                        Column {
                            anchors.fill: parent
                            visible: tabMenu.currentIdentifier === "fluorescence"

                            BoundSlider {
                                id: emissionSlider
                                text: "Emission"
                                minimumValue: 0.01
                                maximumValue: 1000.0
                                target: simulator
                                property: "emissionFactor"
                            }
                        }

                        Column {
                            anchors.fill: parent
                            visible: tabMenu.currentIdentifier === "microscope"

                            BoundSlider {
                                text: "Field of view"
                                minimumValue: 30 / 180
                                maximumValue: 90 / 180
                                target: simulator
                                property: "fieldOfView"
                            }

                            BoundSlider {
                                text: "Lens radius"
                                minimumValue: 0
                                maximumValue: 200
                                target: simulator
                                property: "lensRadius"
                            }

                            BoundSlider {
                                text: "Focal depth"
                                minimumValue: 0
                                maximumValue: 1000
                                target: simulator
                                property: "focalDepth"
                            }
                        }
                    }
                }
            }
        }
    }

    Shortcut {
        sequence: "Ctrl+R"
        onActivated: {
            rebuild()
        }
    }
}
