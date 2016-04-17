#include <QCoreApplication>
#include <QtGui>
#include "sobel.h"
#include "gauss.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QImage qImage("/Users/Dat/Desktop/aaa.jpg");


    Matrix grayscaleMatrix = Matrix::getGrayScaleMatrix(qImage);
    Matrix::exportImage(grayscaleMatrix).save("/Users/Dat/Desktop/myoutputs/gauss_grayscale.png");

    Gauss result = Gauss::Builder(grayscaleMatrix).init().gaussPyramid(grayscaleMatrix).build();

    int count = 0;
    for (const auto& layer : result.getPyramid()) {


        Matrix img = std::get<2>(layer);

        std::string name = "/Users/Dat/Desktop/myoutputs/img_" +
                std::to_string(count++)  +  ".png";

        Matrix::exportImage(img).save(name.c_str());
        printf("exported %s\n",name.c_str());
    }

    return 0;
}

