import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Input 2.0
import Qt3D.Logic 2.0
import QtQml 2.2

Entity {
    id: root

    property Camera camera
    property Entity entity
    property real linearSpeed: 10.0
    property real dragSpeed: 1.0
    property real lookSpeed: 360.0
    property real zoomSpeed: 100.0
    property real zoomLimit: 2.0
    property string mode: "translate"
    property MouseDevice mouseSourceDevice
    property KeyboardDevice keyboardSourceDevice

    components: [
        LogicalDevice {
            actions: [
                Action {
                    id: leftMouseButtonAction
                    ActionInput {
                        sourceDevice: mouseSourceDevice
                        buttons: [MouseEvent.LeftButton]
                    }
                },
                Action {
                    id: rightMouseButtonAction
                    ActionInput {
                        sourceDevice: mouseSourceDevice
                        buttons: [MouseEvent.RightButton]
                    }
                },
                Action {
                    id: middleMouseButtonAction
                    ActionInput {
                        sourceDevice: mouseSourceDevice
                        buttons: [MouseEvent.MiddleButton]
                    }
                },
                Action {
                    id: shiftAction
                    ActionInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Shift]
                    }
                },
                Action {
                    id: sAction
                    onActiveChanged: {
                        if(active) {
                            if(mode === "scale") {
                                mode = "translate"
                            } else {
                                mode = "scale"
                            }
                        }
                    }
                    ActionInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_S]
                    }
                },
                Action {
                    id: rAction
                    onActiveChanged: {
                        if(active) {
                            if(mode === "rotate") {
                                mode = "translate"
                            } else {
                                mode = "rotate"
                            }
                        }
                    }
                    ActionInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_R]
                    }
                },
                Action {
                    id: controlAction
                    ActionInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Control]
                    }
                },
                Action {
                    id: altAction
                    ActionInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Alt]
                    }
                }
            ] // actions

            axes: [
                // Mouse
                Axis {
                    id: mouseXAxis
                    AnalogAxisInput {
                        sourceDevice: mouseSourceDevice
                        axis: MouseDevice.X
                    }
                },
                Axis {
                    id: mouseYAxis
                    AnalogAxisInput {
                        sourceDevice: mouseSourceDevice
                        axis: MouseDevice.Y
                    }
                },

                // Keyboard
                Axis {
                    id: keyboardXAxis
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Left]
                        scale: -1.0
                    }
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Right]
                        scale: 1.0
                    }
                },
                Axis {
                    id: keyboardYAxis
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Up]
                        scale: 1.0
                    }
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Down]
                        scale: -1.0
                    }
                }
            ] // axes
        },
        FrameAction {
            property real timeSinceLastAction: 0.0

            function multiplyQuaternion(q1, q2) {
                return Qt.quaternion(q1.scalar * q2.scalar - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
                                     q1.scalar * q2.x + q1.x * q2.scalar + q1.y * q2.z - q1.z * q2.y,
                                     q1.scalar * q2.y + q1.y * q2.scalar + q1.z * q2.x - q1.x * q2.z,
                                     q1.scalar * q2.z + q1.z * q2.scalar + q1.x * q2.y - q1.y * q2.x);
            }

            function inverseQuaternion(q) {
                return Qt.quaternion(q.scalar, -q.x, -q.y, -q.z)
            }

            onTriggered: {
                if(!root.enabled) {
                    return
                }
                if(!rightMouseButtonAction.active) {
                    timeSinceLastAction += dt
                    return
                }
                if(timeSinceLastAction > 0.1) {
                    timeSinceLastAction = 0
                    return
                }

                if(mode === "translate") {
                    var rightVector = camera.viewVector.crossProduct(camera.upVector).normalized()
                    var upVector = camera.upVector.normalized()
                    var direction = rightVector.times(mouseXAxis.value).plus(upVector.times(mouseYAxis.value))
                    root.entity.transform.translation = root.entity.transform.translation.plus(direction.times(dragSpeed))
                } else if(mode === "scale") {
                    var totalValue = mouseXAxis.value + mouseYAxis.value
                    totalValue *= 0.1
                    root.entity.transform.scale3D = root.entity.transform.scale3D.plus(Qt.vector3d(totalValue, totalValue, totalValue))
                } else if(mode === "rotate") {
                    var totalValue = mouseXAxis.value + mouseYAxis.value
                    totalValue *= 10.0
                    var currentRotation = root.entity.transform.rotation
                    var addedRotation = root.entity.transform.fromAxisAndAngle(camera.viewVector, totalValue)

//                    console.log(addedRotation, inverseQuaternion(addedRotation))

//                    if(currentRotation.x === 0.0 && currentRotation.y === 0.0 && currentRotation.z === 0.0) {
//                        root.entity.transform.rotation = Qt.quaternion(Math.sqrt(2) / 2, 0, 0, Math.sqrt(2) / 2)
//                    } else {
                        root.entity.transform.rotation = multiplyQuaternion(addedRotation, currentRotation)
//                    }
                }
            }
        }
    ] // components
}
