#include "descriptors.h"
#define PI 3.14
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
    printf("init\n");
    printf("descriptor's size: %d - grid's center: %d", this->numberOfBins, this->gridCenter);

    return *this;
}

Descriptors::Builder& Descriptors::Builder::descriptors()
{
    printf("Building Descriptors\n");

    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    Matrix gradientValues = Sobel::Builder().sobelXY(sobelX, sobelY).build().getMatrix();
    Matrix gradientOrientations = Sobel::Builder().gradientOrientiations(sobelX, sobelY).build().getMatrix();

    const int binSize = 2*PI / 8;
    for (const auto& point : this->filteredPoIs) {
        Desciptor descriptor(this->numberOfBins);
        //building descriptor
        for(int i = 0; i < 16; i++){
            for(int j = 0; j < 16; j++){
                int x = std::get<0>(point) - this->gridCenter + i;
                int y = std::get<1>(point) - this->gridCenter + j;
                double gValue = gradientValues.getItensityAt(x, y);
                double gOrientiation = gradientOrientations.getItensityAt(x, y);

                int bin1 = gOrientiation / binSize;
                double bin1Center = bin1 * binSize - binSize / 2;
            }
        }

        this->listOfDesciptors.emplace_back(descriptor);
    }

    return *this;
}


Descriptors Descriptors::Builder::build() const
{
    return Descriptors(this->listOfDesciptors);
}
