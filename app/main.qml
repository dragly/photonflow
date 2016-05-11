import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import Photonflow 1.0
//import SimVis 1.0

ApplicationWindow {
    visible: true
    width: 1280
    height: 1024
    title: qsTr("Photonflow")

    RenderView {
        id: renderView
        anchors.fill: parent
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

    Timer {
        id: timer
        property real lastTime: Date.now()
        property real totalTime: 0
        running: true
        repeat: true
        interval: 16

        onTriggered: {
            console.log("Triggered!")
            totalTime += Date.now() - lastTime
            lastTime = Date.now()
            renderView.integrate()
        }
    }
}

