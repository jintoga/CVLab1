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

   printf("size:%d\n",intensities.capacity());
}

Matrix::Matrix(QImage qImage)
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


int Matrix::getIndex(int row,int col)
{
    return row*width+col;
}

int Matrix::getWidth(){
    return width;
}

int Matrix::getHeight(){
    return height;
}

double Matrix::getItensityAt(int row,int col){
    return intensities[getIndex(row,col)];
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


