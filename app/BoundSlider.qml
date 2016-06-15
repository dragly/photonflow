import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import QtQuick.Layouts 1.1
import Qt.labs.settings 1.0

import Photonflow 1.0

Column {
    id: root
    readonly property bool valid: (root.target && root.property) ? true : false
    readonly property string label: root.text ? root.text + ": " : ""
    property var target
    property string property
    property string text
    property int precision: 2
    property alias minimumValue: slider.minimumValue
    property alias maximumValue: slider.maximumValue
    property alias stepSize: slider.stepSize

    anchors {
        left: parent.left
        right: parent.right
    }

    Text {
        anchors {
            left: parent.left
            right: parent.right
        }

        text: root.valid ? root.label + root.target[root.property].toFixed(precision) : ""
        color: "#ccc"
    }

    Slider {
        id: slider

        anchors {
            left: parent.left
            right: parent.right
        }

        value: root.valid ? root.target[root.property] : 0.0
        stepSize: 0.0
        minimumValue: 0.0
        maximumValue: 1.0

        Binding {
            target: root.valid ? root.target : null
            property: root.valid ? root.property : ""
            value: slider.value
        }
    }

}
