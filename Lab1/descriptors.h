#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>
#include "matrix.h"
#include "pointsofinterest.h"

using Desciptor = std::vector<double>;

class Descriptors
{
public:
    Descriptors();
    class Builder
    {
    public:
        Builder();
        Builder(const Matrix& matrix, const Points& filteredPoIs);
        Builder& init();
    };
};
