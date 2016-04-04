#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>
#include "array.h"


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

            this->pyramid.emplace_back (this->sigma0, this->sigma0, separableFilter(matrix, gauss));
            return *this;
        }

        static std::pair<Array, Array> getGaussSeparable(double sigmaB){
            int kernelSize = int(std::ceil(3 * sigmaB)) * 2 + 1;
            auto separableFilter = std::make_pair<Array, Array>(kernelSize, kernelSize);

            int center = kernelSize / 2;

            for (int i = 0; i < kernelSize; i++) {
                double val =  sqrt(Gauss::gauss(i - center, i - center, sigmaB));
                separableFilter.first.setItem(i,val);
                separableFilter.second.setItem(i,val);
            }
            return separableFilter;

        }

        static Matrix separableFilter(Matrix& matrix, std::pair<Array, Array> filter){
             Matrix tempMat(matrix.getHeight(),matrix.getWidth());
             int center = filter.first.getSize() / 2;
             for (int i = 0; i < matrix.getHeight(); i++) {
                 for (int j = 0; j < matrix.getWidth(); j++) {
                     double result_intensity = 0;
                     for (int y = 0; y < filter.first.getSize(); y++) {
                         auto intensity = matrix.getItensityAt(i, j + y - center);
                         result_intensity += intensity * filter.first.getItem(y);
                     }
                     tempMat.setIntensity(i, j, result_intensity);
                 }
             }

             Matrix result(tempMat.getHeight(), tempMat.getWidth());

             for (int i = 0; i < tempMat.getHeight(); i++) {
                 for (int j = 0; j < tempMat.getWidth(); j++) {
                     double result_intensity = 0;
                     for (int x = 0; x < filter.second.getSize(); x++) {
                         auto intensity = tempMat.getItensityAt(i + x - center, j);
                         result_intensity += intensity * filter.second.getItem(x);
                     }
                     result.setIntensity(i, j, result_intensity);
                 }
             }

             return result;
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
