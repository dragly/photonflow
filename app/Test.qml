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

Entity {
    MouseDevice {
        id: mouseDevice
    }

    components: [
        MouseHandler {
            sourceDevice: mouseDevice
            
            onReleased: {
                console.log(mouse.x, mouse.y)
            }
        }
    ]
}
