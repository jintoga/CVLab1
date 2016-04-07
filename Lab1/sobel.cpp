#include "sobel.h"
#include <QtGui>
#include "matrix.h"

Sobel::Sobel()
{

}

Matrix &Sobel::getMatrix()
{
    printf("matrixsize:%d\n",matrix.getIntensities().capacity());
    return matrix;
}

Matrix Sobel::getGrayScaleMatrix(QImage& qImage) const
{
    Matrix grayScaleMatrix(qImage);
    for (int i = 0; i < grayScaleMatrix.getHeight(); i++) {
        for (int j = 0; j < grayScaleMatrix.getWidth(); j++) {
            auto gray = qGray(qImage.pixel(j,i));
            //printf("%d\n",gray);
            grayScaleMatrix.setIntensity(i,j,intToDouble(gray));

        }
    }


    return grayScaleMatrix;
}

QImage Sobel::exportImage(Matrix& matrix) const
{

    QImage result(matrix.getWidth(),matrix.getHeight(),QImage::Format_RGB32);
    for (int i = 0; i < matrix.getHeight(); i++) {
        for (int j = 0; j < matrix.getWidth(); j++) {
            auto gray = doubleToInt(matrix.getItensityAt(i, j));
            result.setPixel(j, i, qRgb(gray, gray, gray));
        }
    }

    return result;
}

double Sobel::intToDouble(int intensity) const{
    return double(intensity) / 255;
}

int Sobel::doubleToInt(double intensity) const{
    return int(intensity*255);
}

