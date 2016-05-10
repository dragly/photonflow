import SimVis 1.0
import SimVis.ShaderNodes 1.0
import Qt3D.Core 2.0
import Qt3D.Logic 2.0
import Qt3D.Render 2.0
import QtQuick.Scene3D 2.0
import Photonflow 1.0

Scene3D {
    aspects: ["input", "render", "logic"]
    focus: true
    Visualizer {
        id: visualizer
        NeuronSimulator {
            id: simulator
        }
        Cylinders {
            cylinderData: simulator.cylinderData
        }
    }
}
