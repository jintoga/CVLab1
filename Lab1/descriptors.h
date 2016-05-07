#pragma once
#include <QtGui>
#include <cmath>
#include "matrix.h"
#include "sobel.h"
#include "gauss.h"
#include "pointsofinterest.h"

using Desciptor = std::vector<double>;

class Descriptors
{
public:
    Descriptors();
    class Builder
    {
    private:
        Matrix matrix;
        Points filteredPoIs;
    public:
        Builder();
        Builder(const Matrix& matrix, const Points& filteredPoIs);
        Builder& init();
    };
};
