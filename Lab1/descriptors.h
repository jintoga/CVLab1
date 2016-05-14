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
using Orientations = std::vector<double>;
class Descriptors
{
private:
    ListOfDesciptors descriptors;
public:
    Descriptors();
    Descriptors(const ListOfDesciptors& descriptors)
        :descriptors(descriptors)
    {}
    static ResultOfComparision compareDescriptors(const ListOfDesciptors& descriptors1,
                                           const ListOfDesciptors& descriptors2);
    ListOfDesciptors getDescriptors();
    static QImage getMergedMatrix(const Matrix& mat1,
                           const Matrix& mat2,
                           const Points& points1,
                           const Points& points2,
                           const ResultOfComparision& _matches);
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
        Desciptor normalize(const Desciptor& descriptor);
        Builder& descriptors();
        Descriptors build() const;
    };
};
