#pragma once
#include <QtGui>
#include <cmath>
#include "matrix.h"
#include "sobel.h"
#include "gauss.h"
#include "pointsofinterest.h"

using Descriptor = std::vector<double>;
using ListOfDescriptors = std::vector<Descriptor>;
using ResultOfComparision = std::vector<std::pair<int, int>>;
using Orientations = std::vector<double>;
class Descriptors
{
private:
    ListOfDescriptors descriptors;
public:
    Descriptors();
    Descriptors(const ListOfDescriptors& descriptors)
        :descriptors(descriptors)
    {}
    static ResultOfComparision compareDescriptors(const ListOfDescriptors& descriptors1,
                                           const ListOfDescriptors& descriptors2);
    ListOfDescriptors getDescriptors();
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
        ListOfDescriptors listOfDesciptors;
        const int gridCenter = 8;
        const int numberOfBinsPerHistogram = 8;
        const int histogramSize = 4;
        const int numberOfBins = numberOfBinsPerHistogram * histogramSize * histogramSize;
        const int numberOfOrientationBins = 36;
        Matrix gradientValues;
        Matrix gradientOrientations;
    public:
        Builder();
        Builder(const Matrix& matrix, const Points& filteredPoIs);
        Builder& init();
        Descriptor normalize(const Descriptor& descriptor);
        Builder& descriptors();
        Builder& invariantRotationDescriptors();
        Orientations findOrientationBins(const int x, const int y, const double orientationBinSize);
        Descriptor getFinalBins(const Point point, const double mainOrt, const int x, const int y, const double binSize);
        Descriptors build() const;
    };
};
