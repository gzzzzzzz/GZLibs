#pragma once


#include "Runtime/Math/Vector3.h"
#include "Runtime/Math/Quaternion.h"
#include "Runtime/Math/Rect.h"
#include "Runtime/Graphics/Transform.h"
#include "Runtime/Math/Vector2.h"
#include "Runtime/YB/EPSolver.h"
#include "Runtime/Camera/Camera.h"
#include "Runtime/Geometry/Ray.h"
#include "Runtime/Geometry/Plane.h"
#include "Runtime/YB/MathClass.h"


typedef Vector3f Vector3;
typedef Vector2f Vector2;
typedef Quaternionf Quaternion;

void NMSBoxes(const std::vector<Rectf>&bboxes, const std::vector<float>&scores,
    const float score_threshold, const float nms_threshold,
    std::vector<int>& indices, const float eta, const int top_k);
