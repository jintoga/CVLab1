#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>

using Layers = std::tuple<double, double, Matrix>;
using Pyramid = std::vector<Layers>;
using Array = std::vector<double>;
using Filter = std::pair<Array, Array>;
class Gauss
{
private:
    Matrix matrix;
    Pyramid pyramid;
public:
    Gauss();
    Gauss(const Matrix& matrix)
        :matrix(matrix)
    {}
    Gauss(const Pyramid& pyramid)
        :pyramid(pyramid)
    {}
    Pyramid getPyramid() const;
    static double gauss(int x, int y, double sigma);

    class Builder
    {
    private:
        Matrix matrix;
        int s = 3;
        double sigmaA = 0.5;
        double sigma0 = 1.6;
        int n;
        const int minSize = 64;

        Pyramid pyramid;
    public:
        Builder(const Matrix& matrix);
        Builder& init();
        Builder& gaussPyramid(const Matrix& matrix);
        static Filter getGaussSeparable(const double sigmaB);
        static Matrix separableFilter(const Matrix& matrix, const Filter& filter);
        static double getSigmaB(double sigma0, double sigmaA);
        Gauss build() const;
    };
};
