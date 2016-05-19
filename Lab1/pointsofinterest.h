#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>
#include "sobel.h"

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
        int w = 3; //window size
        double threshold = 1;
        int quantity = 300;
        Points pois;
        Points filteredPoIs;
    public:
        Builder();
        Builder(const Matrix& matrix);
        Builder& init();
        Builder& moravec();
        Builder& harris();
        static Matrix opMoravec(const Matrix& matrix, const int w);
        static Matrix opHarris(const Matrix& matrix, const int w);

        static Points findPoI(const Matrix& matrix, const double threshold);
        static Points filterPoI(const Points& points, const unsigned targetQuantity);

        static double getContrast(const Matrix& matrix, const int w,
                                  const int x, const int y,
                                  const int dx, const int dy);
        static double getDistance(const int x1, const int y1,
                                  const int x2, const int y2);
        PointsOfInterest build() const;
    };
};
