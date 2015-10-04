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

    Timer {
        running: true
        repeat: true
        interval: 1000
        onTriggered: {
            renderView.integrate()
        }
    }
}

