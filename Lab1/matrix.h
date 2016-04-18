#pragma once
#include <QtGui>

enum class BorderEffect {
    None,
    Wrap,
    Clamp,
    Mirror
};

class Matrix
{
public:
    Matrix();
    Matrix(const int height,const int width);
    Matrix(const QImage& qImage);
    Matrix(int height,int width,std::vector<double> intensities);
    void setIntensity(int row,int col,double intensity);
    int getIndex(int row,int col) const;

    int getWidth() const;
    int getHeight() const;

    double getItensityAt(int row,int col) const;

    std::vector<double> getIntensities();

    Matrix& normalize();

    static int getRow(int row,int height);

    static int getCol(int col,int width);

    static Matrix getGrayScaleMatrix(const QImage& input);
    static QImage exportImage(const Matrix& matrix);

    static Matrix getDownscaled(const Matrix& mat);
    static Matrix getUpscaled(const Matrix& mat);
private:
    int height;
    int width;
    std::vector<double> intensities;
    BorderEffect borderEffect = BorderEffect::Mirror;

    std::pair<int, int> getClampIndices(int row, int col) const;
    std::pair<int, int> getMirrorIndices(int row, int col) const;
    std::pair<int, int> getWrapIndices(int row, int col) const;
};
