#pragma once
#include <QtGui>
#include <cmath>
#include "matrix.h"
#include "sobel.h"
#include "gauss.h"
#include "pointsofinterest.h"

using Desciptor = std::vector<double>;
using ListOfDesciptors = std::vector<Desciptor>;
using ResultOfComparision = std::vector<std::pair<int, int>>;
class Descriptors
{
private:
    ListOfDesciptors descriptors;
public:
    Descriptors();
    Descriptors(const ListOfDesciptors& descriptors)
        :descriptors(descriptors)
    {}
    ResultOfComparision compareDescriptors(const ListOfDesciptors& descriptors1,
                                           const ListOfDesciptors& descriptors2);
    ListOfDesciptors getDescriptors();
    class Builder
    {
    private:
        Matrix matrix;
        Points filteredPoIs;
        ListOfDesciptors listOfDesciptors;
        const int gridCenter = 8;
        const int numberOfBinsPerHistogram = 8;
        const int histogramSize = 4;
        const int numberOfBins = numberOfBinsPerHistogram * histogramSize * histogramSize;
    public:
        Builder();
        Builder(const Matrix& matrix, const Points& filteredPoIs);
        Builder& init();
        Builder& descriptors();

        Descriptors build() const;
    };
};
