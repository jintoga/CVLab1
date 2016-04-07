#include <QCoreApplication>
#include <QtGui>
#include "sobel.h"
#include "gauss.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QImage qImage("/Users/Dat/Desktop/input.jpg");


    Gauss gauss;
    Matrix grayscaleMatrix = gauss.getGrayScaleMatrix(qImage);
    gauss.exportImage(grayscaleMatrix).save("/Users/Dat/Desktop/myoutputs/gauss_grayscale.png");
    Matrix downscaled = Gauss::getDownscaled(grayscaleMatrix);
    //Matrix upscaled = Gauss::getUpscaled(grayscaleMatrix);

    Gauss result = Gauss::Builder(grayscaleMatrix).init().gaussPyramid(downscaled).build();

    int count = 0;
    for (const auto& layer : result.getPyramid()) {


        Matrix img = std::get<2>(layer);

        std::string name = "/Users/Dat/Desktop/myoutputs/img_" +
                std::to_string(count++)  +  ".png";

        result.exportImage(img).save(name.c_str());
        printf("exported %s\n",name.c_str());
    }

    return 0;
}

