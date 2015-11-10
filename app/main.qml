import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import VSDS 1.0

ApplicationWindow {
    visible: true
    width: 640
    height: 480
    title: qsTr("Hello World")

    RenderView {
        id: renderView
        anchors.fill: parent
    }

    Text {
        color: "white"
        text: Math.round(timer.totalTime, 1)
    }

    MouseArea {
        anchors.fill: parent
        onClicked: timer.running = !timer.running
    }

    Timer {
        id: timer
        property real lastTime: Date.now()
        property real totalTime: 0
        running: true
        repeat: true
        interval: 16

        onTriggered: {
            totalTime += Date.now() - lastTime
            lastTime = Date.now()
            renderView.integrate()
        }
    }
}

