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
using RorationInvariantDescriptor = std::tuple<Descriptor, int, int, double>;
using ListOfRIDescriptors = std::vector<RorationInvariantDescriptor>;
using namespace std;

class Descriptors
{
private:
    ListOfDescriptors descriptors;
    ListOfRIDescriptors riDescriptors;
public:
    Descriptors();
    Descriptors(const ListOfDescriptors& descriptors)
        :descriptors(descriptors)
    {}
    Descriptors(const ListOfRIDescriptors& riDescriptors)
        :riDescriptors(riDescriptors)
    {}
    static ResultOfComparision compareDescriptors(const ListOfDescriptors& descriptors1,
                                           const ListOfDescriptors& descriptors2);
    ListOfDescriptors getDescriptors(boolean isRotationInvariant);
    ListOfRIDescriptors getRIDescriptors();

    static QImage getMergedMatrix(const Matrix& mat1,
                           const Matrix& mat2,
                           const Points& points1,
                           const Points& points2,
                           const ResultOfComparision& _matches);
    static QImage getMergedMatrix(const Matrix& mat1,
                           const Matrix& mat2,
                           const ListOfRIDescriptors& d1,
                           const ListOfRIDescriptors& d2,
                           const ResultOfComparision& _matches);
    class Builder
    {
    private:
        Matrix matrix;
        Points filteredPoIs; 
        ListOfRIDescriptors listOfRIDescriptors;
        const int gridCenter = 8;
        const int binsPerHistogram = 8;
        const int binsOfWideHistogram = 36;
        const int histogramSize = 4;
        const int numberOfBins = binsPerHistogram * histogramSize * histogramSize;
    public:
        Builder();
        Builder(const Matrix& matrix, const Points& filteredPoIs);
        Builder& init();
        Descriptor normalize(const Descriptor& descriptor); 
        double interpolation(const double x2, const double y1, const double y2, const double y3);
        Builder& descriptors();
        Descriptors build() const;
    };
};
