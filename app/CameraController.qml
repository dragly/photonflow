import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Input 2.0
import Qt3D.Logic 2.0
import QtQml 2.2

Entity {
    id: root
    property Camera camera
    property real linearSpeed: 10.0
    property real lookSpeed: 100.0
    property real zoomSpeed: 20.0
    property real zoomLimit: 2.0
    property MouseDevice mouseSourceDevice
    property KeyboardDevice keyboardSourceDevice

    QtObject {
        id: d
        readonly property vector3d firstPersonUp: Qt.vector3d(0, 1, 0)
    }

    function zoomDistance(firstPoint, secondPoint) {
        var u = secondPoint.minus(firstPoint); u = u.times(u);
        return u.x + u.y + u.z;
    }

    components: [
        LogicalDevice {
            actions: [
                Action {
                    id: noMouseButtonAction
                    ActionInput {
                        sourceDevice: mouseSourceDevice
                        buttons: [MouseEvent.NoButton]
                    }
                },
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
                    id: keyboardZAxis
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Up]
                        scale: d.shiftPressed ? 1.0 : 0.0
                    }
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Down]
                        scale: d.shiftPressed ? -1.0 : 0.0
                    }
                },
                Axis {
                    id: keyboardYAxis
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Up]
                        scale: d.shiftPressed ? 0.0 : 1.0
                    }
                    ButtonAxisInput {
                        sourceDevice: keyboardSourceDevice
                        buttons: [Qt.Key_Down]
                        scale: d.shiftPressed ? 0.0 : -1.0
                    }
                }
            ] // axes
        },

        FrameAction {
            onTriggered: {
                if(!root.enabled) {
                    return
                }

                // The time difference since the last frame is passed in as the
                // argument dt. It is a floating point value in units of seconds.
                if (leftMouseButtonAction.active) {
                    root.camera.panAboutViewCenter(-mouseXAxis.value * lookSpeed, d.firstPersonUp);
                    root.camera.tiltAboutViewCenter(-mouseYAxis.value * lookSpeed);
                } else if(middleMouseButtonAction.active) {
                    var fov = root.camera.fieldOfView
                    fov += mouseYAxis.value * zoomSpeed
                    fov = Math.max(10.0, Math.min(160.0, fov))
                    root.camera.fieldOfView = fov
                }
                root.camera.translate(Qt.vector3d(keyboardXAxis.value * 0.3 * linearSpeed * dt, 0.0, keyboardYAxis.value * linearSpeed * dt))
            }
        }
    ] // components
}
