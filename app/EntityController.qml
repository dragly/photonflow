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
    property real dragSpeed: 10.0
    property real lookSpeed: 360.0
    property real zoomSpeed: 100.0
    property real zoomLimit: 2.0

    onEnabledChanged: {
        console.log("Controller enabled", enabled)
    }

    KeyboardDevice {
        id: keyboardSourceDevice
    }

    MouseDevice {
        id: mouseSourceDevice
        sensitivity: 0.1
    }

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
            property bool buttonWasDownPreviousFrame: false
            onTriggered: {
                if(!root.enabled) {
                    return
                }
                if(!rightMouseButtonAction.active) {
                    return
                }

                // The time difference since the last frame is passed in as the
                // argument dt. It is a floating point value in units of seconds.
                var rightVector = camera.viewVector.crossProduct(camera.upVector).normalized()
                var upVector = camera.upVector.normalized()
                var direction = rightVector.times(mouseXAxis.value).plus(upVector.times(mouseYAxis.value))
                root.entity.transform.translation = root.entity.transform.translation.plus(direction.times(dragSpeed * dt))
            }
        }
    ] // components
}
