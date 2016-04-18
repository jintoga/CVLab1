#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>

using Point = std::tuple<int, int, double>;
using Points = std::vector<Point>;

class PointsOfInterest
{
private:
    Matrix matrix;
    Points pois;
    Points filteredPoIs;
public:
    PointsOfInterest();
    PointsOfInterest(const Matrix& matrix,
                     const Points& pois,
                     const Points& filteredPoIs)
        :matrix(matrix)
        ,pois(pois)
        ,filteredPoIs(filteredPoIs)
    {}
    Matrix getMatrix();
    Points getPoIs();
    Points getFilteredPoIs();

    static QImage markPoints(const Matrix& matrix, const Points& points);
    class Builder
    {
    private:
        Matrix matrix;
        int w = 3;
        double threshold = 0.15;
        int p = 5;
        int quantity = 100;
        Points pois;
        Points filteredPoIs;
    public:
        Builder();
        Builder(const Matrix& matrix);
        Builder& init();
        Builder& moravec();
        static Matrix opMoravec(const Matrix& matrix, const int w);

        static Points findPoI(const Matrix& matrix, const double threshold, const int p);
        static Points filterPoI(const Points& points, const unsigned targetQuantity);

        static double getC(const Matrix& matrix, const int w,
                           const int x, const int y,
                           const int dx, const int dy);
        static double getDistance(const int x1, const int y1,
                                  const int x2, const int y2);
        PointsOfInterest build() const;
    };
};