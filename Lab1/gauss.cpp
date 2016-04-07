#include "gauss.h"

Gauss::Gauss()
{

}

Pyramid Gauss::getPyramid() const
{
    return pyramid;
}

Matrix Gauss::getGrayScaleMatrix(QImage& qImage) const
{
    Matrix grayScaleMatrix(qImage);
    for (int i = 0; i < grayScaleMatrix.getHeight(); i++) {
        for (int j = 0; j < grayScaleMatrix.getWidth(); j++) {
            auto gray = qGray(qImage.pixel(j,i));
            grayScaleMatrix.setIntensity(i,j,intToDouble(gray));

        }
    }
    return grayScaleMatrix;
}

Matrix Gauss::getDownscaled(Matrix& mat)
{
    Matrix result(mat.getHeight()/2, mat.getWidth()/2);

    for(int i = 0; i < result.getHeight();i++){
        for(int j = 0;j < result.getWidth();j++){
            auto gradient = mat.getItensityAt(i * 2, j * 2);
            result.setIntensity(i, j, gradient);
        }
    }

    printf("downscaled\n");
    return result;
}

Matrix Gauss::getUpscaled(Matrix& mat)
{
    Matrix result(mat.getHeight()*2, mat.getWidth()*2);

    for(int i = 0; i < mat.getHeight();i++){
        for(int j = 0;j < mat.getWidth();j++){
            result.setIntensity(i * 2, j * 2, mat.getItensityAt(i, j));
            result.setIntensity(i * 2, j * 2 + 1, mat.getItensityAt(double(i), j + 0.5));
            result.setIntensity(i * 2 + 1, j * 2, mat.getItensityAt(i + 0.5, double(j)));
            result.setIntensity(i * 2 + 1, j * 2 + 1, mat.getItensityAt(i + 0.5, j + 0.5));
        }
    }

    printf("upscaled\n");
    return result;
}


double Gauss::gauss(int x, int y, double sigma){
    double val = 2 * sigma * sigma;
    return 1.0 / (3.14 * val) * exp(-(x * x + y * y) / val);
}

QImage Gauss::exportImage(Matrix& matrix) const
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

double Gauss::intToDouble(int intensity) const{
    return double(intensity) / 255;
}

int Gauss::doubleToInt(double intensity) const{
    return int(intensity*255);
}
