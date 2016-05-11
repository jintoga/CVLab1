#include <QCoreApplication>
#include <QtGui>
#include "sobel.h"
#include "gauss.h"
#include "matrix.h"
#include "pointsofinterest.h"
#include "descriptors.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QImage qImage("/Users/Dat/Desktop/lena.jpg");


    Matrix grayscaleMatrix = Matrix::getGrayScaleMatrix(qImage);
    Matrix::exportImage(grayscaleMatrix).save("/Users/Dat/Desktop/myoutputs/grayscale.png");


    PointsOfInterest harris = PointsOfInterest::Builder(grayscaleMatrix).init().harris().build();
    printf("points of interest filtered: %d\n\n", harris.getFilteredPoIs().size());
    Matrix::exportImage(harris.getMatrix().normalize()).save("/Users/Dat/Desktop/myoutputs/harris.png");
    harris.markPoints(grayscaleMatrix, harris.getFilteredPoIs()).save("/Users/Dat/Desktop/myoutputs/harris_filtered_pois.png");

    Descriptors descriptors = Descriptors::Builder(grayscaleMatrix, harris.getFilteredPoIs()).init().descriptors().build();
    printf("descriptors: %d\n", descriptors.getDescriptors().size());
    Matrix::exportImage(harris.getMatrix().normalize()).save("/Users/Dat/Desktop/myoutputs/harris.png");

    return 0;
}

