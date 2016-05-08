#pragma once
#include <QtGui>
#include <cmath>
#include "matrix.h"
#include "sobel.h"
#include "gauss.h"
#include "pointsofinterest.h"

using Desciptor = std::vector<double>;
using ListOfDesciptors = std::vector<Desciptor>;

class Descriptors
{
private:
    ListOfDesciptors descriptors;
public:
    Descriptors();
    Descriptors(const ListOfDesciptors& descriptors)
        :descriptors(descriptors)
    {}
    ListOfDesciptors getDescriptors();
    class Builder
    {
    private:
        Matrix matrix;
        Points filteredPoIs;
        ListOfDesciptors listOfDesciptors;
        const int numberOfBins = 8 * 4 * 4;
        const int gridCenter = 8;
    public:
        Builder();
        Builder(const Matrix& matrix, const Points& filteredPoIs);
        Builder& init();
        Builder& descriptors();

        Descriptors build() const;
    };
};
