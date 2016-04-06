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
    Gauss(std::vector<std::tuple<double, double, Matrix>>& pyramid)
        :pyramid(pyramid)
    {};
    Matrix getGrayScaleMatrix(QImage& input);
    double intToDouble(int intensity) const;
    int doubleToInt(double intensity) const;
    QImage exportImage(Matrix& matrix) const;
    std::vector<std::tuple<double, double, Matrix>> getPyramid();
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
            this->n = this->n + 1;
            printf("const: s = %d, sigmaA = %lf, sigma0 = %lf, n = %d\n", this->s, this->sigmaA, this->sigma0, n);
            return *this;
        }

        Builder& gaussPyramid(Matrix& matrix){
            auto sigmaB = getSigmaB(this->sigma0, this->sigmaA);
            auto gauss = getGaussSeparable(sigmaB);

            this->pyramid.emplace_back (this->sigma0, this->sigma0, separableFilter(matrix, gauss));

            auto k = pow(2.0, 1.0 / this->s);

            std::vector<std::pair<Array, Array>> filters;
            auto oldSigma = this->sigma0;
            for (int i = 0; i < this->s; i++) {
                auto newSigma = oldSigma * k;
                auto sigma = getSigmaB(newSigma, oldSigma);
                filters.push_back(getGaussSeparable(sigma));
                oldSigma = newSigma;
            }

            for (int oct = 0; oct < this->n; oct++) {
                 for (int i = 0; i < this->s; i++) {
                    std::tuple<double, double, Matrix>& layer = this->pyramid.back();
                    double s0 = std::get<0>(layer);
                    double s1 = std::get<1>(layer);

                    Matrix& image = std::get<2>(layer);

                    this->pyramid.emplace_back(s0*k, s1*k, separableFilter(image, filters[size_t(i)]));
                 }

                 if (oct != this->n - 1) {

                     std::tuple<double, double, Matrix>& layer = this->pyramid.back();
                     auto s = std::get<1>(layer);
                     Matrix& image = std::get<2>(layer);

                     this->pyramid.emplace_back(this->sigma0, s, getDownscale(image));
                 }
            }
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

        static Matrix separableFilter(Matrix& matrix,const std::pair<Array, Array> filter){
             Matrix tempMat(matrix.getHeight(),matrix.getWidth());
             int center = filter.first.getSize() / 2;
             for (int i = 0; i < matrix.getHeight(); i++) {
                 for (int j = 0; j < matrix.getWidth(); j++) {
                     double result_intensity = 0;
                     for (int y = 0; y < filter.first.getSize(); y++) {
                         int jj = j + y - center;
                         jj = Matrix::getCol(jj,tempMat.getWidth());
                         auto intensity = matrix.getItensityAt(i, jj);
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
                         int ii = i + x - center;
                         ii = Matrix::getRow(ii,tempMat.getHeight());
                         auto intensity = tempMat.getItensityAt(ii, j);
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
            return Gauss(this->pyramid);
        }
    };
};
