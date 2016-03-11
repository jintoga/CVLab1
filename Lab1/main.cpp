#include <QCoreApplication>
#include <QtGui>
#include "sobel.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QImage qImage("/Users/Dat/Desktop/input.jpg");

    Sobel sobel;
    Matrix grayscaleMatrix = sobel.getGrayScaleMatrix(qImage);
    sobel.exportImage(grayscaleMatrix).save("/Users/Dat/Desktop/outputs/grayscale.png");


    Matrix matSobelX = Sobel::Builder(grayscaleMatrix).sobelX().build().getMatrix();
    Matrix matSobelY = Sobel::Builder(grayscaleMatrix).sobelY().build().getMatrix();

    sobel.exportImage(matSobelX.normalize()).save("/Users/Dat/Desktop/outputs/sobel_x.png");
    sobel.exportImage(matSobelY.normalize()).save("/Users/Dat/Desktop/outputs/sobel_y.png");

    Matrix x = Sobel::Builder(grayscaleMatrix).sobelX().build().getMatrix();
    Matrix y = Sobel::Builder(grayscaleMatrix).sobelY().build().getMatrix();

    Matrix matSobelXY = Sobel::Builder(grayscaleMatrix).sobelXY(x,y).build().getMatrix();
    sobel.exportImage(matSobelXY.normalize()).save("/Users/Dat/Desktop/outputs/sobel_xy.png");


    return 0;
}

