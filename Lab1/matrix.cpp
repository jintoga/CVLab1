#include "matrix.h"
#include <cassert>

Matrix::Matrix()
{

}

Matrix::Matrix(int height,int width)
    :height(height)
    ,width(width)
    ,intensities(height * width)
{
}

Matrix::Matrix(const QImage& qImage)
    :Matrix(qImage.height(),qImage.width())
{
}

Matrix::Matrix(int height,int width,std::vector<double> intensities)
    :Matrix(height,width)
{
    this->intensities = intensities;
}

void Matrix::setIntensity(int row,int col,double val)
{
    intensities[getIndex(row,col)] = val;
}


int Matrix::getIndex(int row,int col) const
{
    return row*width+col;
}

int Matrix::getWidth() const{
    return width;
}

int Matrix::getHeight() const{
    return height;
}

double Matrix::getItensityAt(int row,int col) const{
    return intensities[getIndex(row,col)];
}

int Matrix::getRow(int row,int height){
    int res_row;
    if (row < 0) {
        res_row = 0;
    } else if (row >= height) {
        res_row = height - 1;
    } else {
        res_row = row;
    }
    return res_row;
}

int Matrix::getCol(int col,int width){
    int res_col;
    if (col < 0) {
        res_col = 0;
    } else if (col >= width) {
        res_col = width - 1;
    } else {
        res_col = col;
    }
    return res_col;
}

std::vector<double> Matrix::getIntensities(){
    return intensities;
}

Matrix& Matrix::normalize(){
    auto minmax = std::minmax_element(intensities.begin(), intensities.end());
    auto min_intensity = *minmax.first;
    auto max_intensity = *minmax.second;
    std::transform(intensities.begin(), intensities.end(), intensities.begin(),
    [=](const auto& intensity)
    {
        return (intensity - min_intensity) / double(max_intensity - min_intensity);

    });
    return *this;
}

double intToDouble(const int intensity){
    return double(intensity) / 255;
}

int doubleToInt(const double intensity){
    return int(intensity*255);
}


Matrix Matrix::getGrayScaleMatrix(const QImage& qImage)
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


QImage Matrix::exportImage(const Matrix& matrix)
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


Matrix Matrix::getDownscaled(const Matrix& mat)
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

Matrix Matrix::getUpscaled(const Matrix& mat)
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


