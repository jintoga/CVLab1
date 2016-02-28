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

    sobel.exportImage(grayscaleMatrix).save("/Users/Dat/Desktop/outputs/grayscale.jpg");

    Sobel::Builder builderSobelX(grayscaleMatrix);
    Sobel::Builder builderSobelY(grayscaleMatrix);

    Matrix matSobelX = builderSobelX.sobelX().build().getMatrix();
    Matrix matSobelY = builderSobelY.sobelY().build().getMatrix();

    sobel.exportImage(matSobelX.normalize()).save("/Users/Dat/Desktop/outputs/sobel_x.jpg");
    sobel.exportImage(matSobelY.normalize()).save("/Users/Dat/Desktop/outputs/sobel_y.jpg");

    Sobel::Builder builderSobelXY(grayscaleMatrix);

    Matrix x = builderSobelX.sobelX().build().getMatrix();
    Matrix y = builderSobelY.sobelY().build().getMatrix();

    Matrix matSobelXY = builderSobelXY.sobelXY(x,y).build().getMatrix();
    sobel.exportImage(matSobelXY.normalize()).save("/Users/Dat/Desktop/outputs/sobel_xy.jpg");

    Sobel::Builder builderSobel(grayscaleMatrix);

    Matrix matSobel = builderSobel.sobelX().sobelY().build().getMatrix();

    sobel.exportImage(matSobel.normalize()).save("/Users/Dat/Desktop/outputs/sobel.jpg");

    return 0;
}

