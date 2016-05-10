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
    printf("descriptor's size: %d - grid's center: %d\n", this->numberOfBins, this->gridCenter);

    return *this;
}

Descriptors::Builder& Descriptors::Builder::descriptors()
{
    printf("Building Descriptors\n");

    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    Matrix gradientValues = Sobel::Builder().sobelXY(sobelX, sobelY).build().getMatrix();
    Matrix gradientOrientations = Sobel::Builder().gradientOrientiations(sobelX, sobelY).build().getMatrix();

    const double binSize = 3.14*2 / this->numberOfBinsPerHistogram;

    for (const auto& point : this->filteredPoIs) {
        Desciptor descriptor(this->numberOfBins);
        int x = std::get<0>(point) - this->gridCenter;
        int y = std::get<1>(point) - this->gridCenter;

        //building descriptor's histograms
        for(int i = 0; i < 16; i++){
            for(int j = 0; j < 16; j++){


                double gValue = gradientValues.getItensityAt(x + i, y + j);
                double gOrientation = gradientOrientations.getItensityAt(x + i, y + j);

                //indexing bins
                int bin1Index = gOrientation / binSize;
                double bin1Center = bin1Index * binSize + binSize / 2;

                //int bin2Index = gOrientation >= bin1Center ? bin1Index + 1 : bin1Index - 1;
                int bin2Index = bin1Index + 1;
                if(gOrientation < bin1Center)
                    bin2Index  = bin1Index - 1;
                //check for histogram's edges
                bin2Index = (bin2Index + this->numberOfBinsPerHistogram) % this->numberOfBinsPerHistogram;

                //calculating distance to center
                double bin1Dist = abs(gOrientation - bin1Center);
                double bin2Dist = binSize - bin1Dist;

                //get current histogram's index
                int curHistIndex = (i / this->histogramSize * 4 + j / this->histogramSize) * 4;

                descriptor[curHistIndex + bin1Index] += gValue * (1 - bin1Dist / binSize) ;
                descriptor[curHistIndex + bin2Index] += gValue * (1 - bin2Dist / binSize) ;

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
