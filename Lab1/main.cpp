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

    QImage qImage1("/Users/Dat/Desktop/lena1.jpg");

    QImage qImage2("/Users/Dat/Desktop/lena2.jpg");

    Matrix grayscaleMatrix1 = Matrix::getGrayScaleMatrix(qImage1);
    Matrix::exportImage(grayscaleMatrix1).save("/Users/Dat/Desktop/myoutputs/grayscale1.png");

    Matrix grayscaleMatrix2 = Matrix::getGrayScaleMatrix(qImage2);
    Matrix::exportImage(grayscaleMatrix2).save("/Users/Dat/Desktop/myoutputs/grayscale2.png");


    PointsOfInterest harris = PointsOfInterest::Builder(grayscaleMatrix1).init().harris().build();
    printf("points of interest filtered: %d\n\n", harris.getFilteredPoIs().size());
    Matrix::exportImage(harris.getMatrix().normalize()).save("/Users/Dat/Desktop/myoutputs/harris1.png");
    harris.markPoints(grayscaleMatrix1, harris.getFilteredPoIs()).save("/Users/Dat/Desktop/myoutputs/harris_filtered_pois1.png");

    Descriptors descriptors = Descriptors::Builder(grayscaleMatrix1, harris.getFilteredPoIs()).init().descriptors().build();
    printf("descriptors: %d\n", descriptors.getDescriptors().size());

    PointsOfInterest harris2 = PointsOfInterest::Builder(grayscaleMatrix2).init().harris().build();
    printf("points of interest filtered: %d\n\n", harris2.getFilteredPoIs().size());
    Matrix::exportImage(harris2.getMatrix().normalize()).save("/Users/Dat/Desktop/myoutputs/harris2.png");
    harris.markPoints(grayscaleMatrix2, harris2.getFilteredPoIs()).save("/Users/Dat/Desktop/myoutputs/harris_filtered_pois2.png");

    Descriptors descriptors2 = Descriptors::Builder(grayscaleMatrix2, harris2.getFilteredPoIs()).init().descriptors().build();
    printf("descriptors: %d\n", descriptors2.getDescriptors().size());

    Descriptors d;
    ResultOfComparision matches = d.compareDescriptors(descriptors.getDescriptors(), descriptors2.getDescriptors());

    d.getMergedMatrix(grayscaleMatrix1,
                      grayscaleMatrix2,
                      harris.getFilteredPoIs(),
                      harris2.getFilteredPoIs(),
                      matches).save("/Users/Dat/Desktop/myoutputs/merged.png");
    return 0;
}

