#pragma once
#include <QtGui>
#include "matrix.h"
#include <cmath>



class Sobel
{

private:
    Matrix matrix;

public:
    Sobel();
    Sobel(Matrix& matrix)
        :matrix(matrix)
    {}
    Matrix& getMatrix();

    class Builder
    {
    private:

        const int size = 3;

        const std::vector<double> sobel_x_vector{
            -1, 0, 1,
            -2, 0, 2,
            -1, 0, 1
        };
        const std::vector<double> sobel_y_vector{
            -1, -2, -1,
             0,  0,  0,
             1,  2,  1
        };
    public:
        Builder()
        {}
        Builder(const Matrix& matrix)
            :matrix(matrix)
        {}
        Builder& sobelX()
        {
            printf("sobelX\n");
            Matrix sob(size,size,sobel_x_vector);

            this->matrix = convolution(this->matrix,sob);

            return *this;
        }

        Builder& sobelY()
        {
            printf("sobelY\n");
            Matrix sob(size,size,sobel_y_vector);

            this->matrix = convolution(this->matrix,sob);

            return *this;
        }

        Builder& sobelXY(const Matrix& mat1, const Matrix& mat2)
        {
            printf("sobelXY\n");
            Matrix result(mat1.getHeight(),mat1.getWidth());

            for(int i = 0; i < result.getHeight();i++){
                for(int j = 0;j < result.getWidth();j++){
                    auto gradient = getGradient(mat1.getItensityAt(i, j), mat2.getItensityAt(i, j));
                    result.setIntensity(i, j, gradient);
                }
            }

            this->matrix = result;

            return *this;
        }

        Builder& gradientOrientiations(const Matrix& mat1, const Matrix& mat2)
        {
            printf("gradientDirections\n");
            Matrix result(mat1.getHeight(),mat1.getWidth());

            for(int i = 0; i < result.getHeight();i++){
                for(int j = 0;j < result.getWidth();j++){
                    auto gradient = getGradientDirection(mat1.getItensityAt(i, j), mat2.getItensityAt(i, j));
                    result.setIntensity(i, j, gradient);
                }
            }

            this->matrix = result;

            return *this;
        }

        Sobel build()
        {
            return Sobel(this->matrix);
        }

        static Matrix convolution(const Matrix& mat1, const Matrix& mat2)
        {
            Matrix result(mat1.getHeight(),mat1.getWidth());


            int centerX = mat2.getHeight() / 2;
            int centerY = mat2.getWidth() / 2;

            for(int i=0; i < mat1.getHeight(); i++){
                for(int j = 0; j < mat1.getWidth(); j++){
                    double result_intensity = 0;
                    for(int x = 0; x < mat2.getHeight(); x++){
                        for(int y = 0; y < mat2.getWidth(); y++){
                            int ii = i + x - centerX;
                            int jj = j + y - centerY;
                            if( ii >= 0 && ii < mat1.getHeight() && jj >= 0 && jj <  mat1.getWidth())
                                result_intensity += mat2.getItensityAt(x, y) * mat1.getItensityAt(ii,jj);
                            //edge effect
                            else{
                                ii = Matrix::getRow(ii,mat1.getHeight());
                                jj = Matrix::getCol(jj,mat1.getWidth());
                                result_intensity += mat2.getItensityAt(x, y) * mat1.getItensityAt(ii,jj);
                            }
                        }
                    }
                    result.setIntensity(i,j,result_intensity);
                }
            }

            return result;
        }

        double getGradient(double x, double y){
            return sqrt(x * x + y * y);
        }

        double getGradientDirection(double x, double y){
            double res = atan2(y, x);
            return res >= 0 ? res : res + M_PI * 2;
        }

    private:
        Matrix matrix;

    };


};



