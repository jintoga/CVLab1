#pragma once
#include <QtGui>

class Matrix
{
public:
    Matrix();
    Matrix(int height,int width);
    Matrix(QImage qImage);
    Matrix(int height,int width,std::vector<double> intensities);
    void setIntensity(int row,int col,double intensity);
    int getIndex(int row,int col);

    int getWidth();
    int getHeight();

    double getItensityAt(int row,int col);

    std::vector<double> getIntensities();

    Matrix& normalize();

private:
    int height;
    int width;
    std::vector<double> intensities;

};
