import QtQuick 2.5
import QtQuick.Controls 1.4
import QtQuick.Dialogs 1.2
import QtQuick.Layouts 1.1
import QtGraphicalEffects 1.0
import Qt.labs.settings 1.0

import Photonflow 1.0

Item {
    property bool toggled
    property alias source: image.source
    
    Image {
        id: image
        anchors.fill: parent
        fillMode: Image.PreserveAspectFit
    }
    ColorOverlay {
        anchors.fill: parent
        source: image
        color: "#888"
        visible: !toggled
    }
    Glow {
        anchors.fill: parent
        visible: toggled
        source: image
        samples: 12
        color: "#3377ef"
    }
}
