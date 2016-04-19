#include <QCoreApplication>
#include <QtGui>
#include "sobel.h"
#include "gauss.h"
#include "matrix.h"
#include "pointsofinterest.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QImage qImage("/Users/Dat/Desktop/lena_std.tif");


    Matrix grayscaleMatrix = Matrix::getGrayScaleMatrix(qImage);
    Matrix::exportImage(grayscaleMatrix).save("/Users/Dat/Desktop/myoutputs/grayscale.png");

    PointsOfInterest moravec = PointsOfInterest::Builder(grayscaleMatrix).init().moravec().build();
    printf("points of interest: %d\n", moravec.getPoIs().size());
    printf("points of interest filtered: %d\n\n", moravec.getFilteredPoIs().size());
    Matrix::exportImage(moravec.getMatrix()).save("/Users/Dat/Desktop/myoutputs/moravec.png");
    moravec.markPoints(grayscaleMatrix, moravec.getPoIs()).save("/Users/Dat/Desktop/myoutputs/moravec_pois.png");
    moravec.markPoints(grayscaleMatrix, moravec.getFilteredPoIs()).save("/Users/Dat/Desktop/myoutputs/moravec_filtered_pois.png");

    PointsOfInterest harris = PointsOfInterest::Builder(grayscaleMatrix).init().harris().build();
    printf("points of interest: %d\n", harris.getPoIs().size());
    printf("points of interest filtered: %d\n\n", harris.getFilteredPoIs().size());
    Matrix::exportImage(harris.getMatrix().normalize()).save("/Users/Dat/Desktop/myoutputs/harris.png");
    harris.markPoints(grayscaleMatrix, harris.getPoIs()).save("/Users/Dat/Desktop/myoutputs/harris_pois.png");
    harris.markPoints(grayscaleMatrix, harris.getFilteredPoIs()).save("/Users/Dat/Desktop/myoutputs/harris_filtered_pois.png");


    return 0;
}

