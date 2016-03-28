#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>


class Gauss
{
private:
    Matrix matrix;

    std::tuple<double, double, Matrix> pyramidLayers;
    std::vector<std::tuple<double, double, Matrix>> pyramid;
public:
    Gauss();
    Gauss(Matrix& matrix)
        :matrix(matrix)
    {};

    Matrix getGrayScaleMatrix(QImage& input);
    double intToDouble(int intensity) const;
    int doubleToInt(double intensity) const;
    QImage exportImage(Matrix& matrix) const;
    Matrix& getMatrix();
    static Matrix getDownscale(Matrix& mat);
    static double gauss(int x, int y, double sigma);

    class Builder
    {
    private:
        Matrix matrix;
        int s = 3;
        double sigmaA = 0.5;
        double sigma0 = 1.6;
        int n;
        const int min_size = 64;

        std::vector<std::tuple<double, double, Matrix>> pyramid;
    public:
        Builder(Matrix& matrix)
            :matrix(matrix)
        {}

        Builder& init()
        {
            printf("init\n");
            this->n = int(std::min(std::log2(double(matrix.getHeight()) / min_size),
                              std::log2(double(matrix.getWidth()) / min_size))) + 1;
            printf("const: s = %d, sigmaA = %lf, sigma0 = %lf, n = %d\n", this->s, this->sigmaA, this->sigma0, n);
            return *this;
        }

        Builder& gaussPyramid(Matrix& matrix){
            auto sigmaB = getSigmaB(this->sigma0, this->sigmaA);
            auto gauss = getGaussSeparable(sigmaB);
            return *this;
        }

        static std::vector<std::pair<int, int>> getGaussSeparable(double sigmaB){
            int kernelSize = int(std::ceil(3 * sigmaB)) * 2 + 1;
            std::vector<std::pair<int, int>> separableFilter;

            int half = kernelSize / 2;

            for (int i = 0; i < kernelSize; i++) {
                separableFilter[i].first = separableFilter[i].second = sqrt(Gauss::gauss(i - half, i - half, sigmaB));
            }
            return separableFilter;

        }

        static double getSigmaB(double sigma0, double sigmaA) {
            return sqrt(sigma0 * sigma0 - sigmaA * sigmaA);
        }

        Gauss build()
        {
            return Gauss(this->matrix);
        }
    };
};
