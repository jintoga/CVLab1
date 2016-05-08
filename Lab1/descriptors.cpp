#include "descriptors.h"

Descriptors::Descriptors()
{

}

ListOfDesciptors Descriptors::getDescriptors()
{
    return this->descriptors;
}

Descriptors::Builder::Builder()
{}

Descriptors::Builder::Builder(const Matrix& matrix, const Points& filteredPoIs)
    :matrix(matrix)
    ,filteredPoIs(filteredPoIs)
{}

Descriptors::Builder& Descriptors::Builder::init()
{

    return *this;
}

Descriptors::Builder& Descriptors::Builder::descriptors()
{
    printf("Descriptors\n");

    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    Matrix gradientValues = Sobel::Builder(matrix).sobelXY(sobelX, sobelY).build().getMatrix();
    Matrix gradientDirections = Sobel::Builder(matrix).gradientDirections(sobelX, sobelY).build().getMatrix();

    for (const auto& point : this->filteredPoIs) {

    }

    return *this;
}


Descriptors Descriptors::Builder::build() const
{
    return Descriptors(this->listOfDesciptors);
}
