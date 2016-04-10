#include "sobel.h"
#include <QtGui>
#include "matrix.h"

Sobel::Sobel()
{

}

Matrix &Sobel::getMatrix()
{
    printf("matrixsize:%d\n",matrix.getIntensities().capacity());
    return matrix;
}


