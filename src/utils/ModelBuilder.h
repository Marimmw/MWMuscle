#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "utils/MWMath.h"

// THIS FILE CONTAINS HELPER FUNCTIONS FOR CREATING MODELS

struct JoingAngleDef
{
    int steps;
    double finalAngleDeg;
};

inline void createJointMovementVector(std::vector<double>& jointAngles, const std::vector<JoingAngleDef>& jointAngleDefs = {{10, 90.0}}) {
    double currentAngle;
    if (jointAngles.empty()){
        currentAngle = 0.0;
    }
    else{
        currentAngle = jointAngles.back();
    }
    for (const JoingAngleDef& def : jointAngleDefs) {
        double stepSize = (def.finalAngleDeg - currentAngle) / def.steps;
        for (int i = 0; i < def.steps; ++i) {
            currentAngle += stepSize;
            jointAngles.push_back(currentAngle);
        }
    }
}