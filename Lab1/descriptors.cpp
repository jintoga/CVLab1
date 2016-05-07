#include "descriptors.h"

Descriptors::Descriptors()
{

}

Descriptors::Builder::Builder()
{}

Descriptors::Builder::Builder(const Matrix& matrix, const Points& filteredPoIs)
    :matrix(matrix)
    ,filteredPoIs(filteredPoIs)
{}

Descriptors::Builder& Descriptors::Builder::init()
{

    Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    return *this;
}

