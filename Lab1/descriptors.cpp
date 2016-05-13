#include "descriptors.h"

Descriptors::Descriptors()
{

}

ListOfDesciptors Descriptors::getDescriptors()
{
    return this->descriptors;
}

static double getDistance(const Desciptor& desc1, const Desciptor& desc2) {
    double sum = 0;
    for (unsigned i = 0; i < desc1.size(); i++) {
        sum += (desc1[i] - desc2[i])*(desc1[i] - desc2[i]);
    }
    return sqrt(sum);
}

static unsigned getIndexOfNearest(const Desciptor& desc, const ListOfDesciptors& descriptors) {

    int minIndex = 0;
    auto minDistance = getDistance(desc, descriptors[minIndex]);

    for (unsigned i = 1; i < descriptors.size(); i++) {
        auto distance = getDistance(desc, descriptors[i]);
        if (distance < minDistance) {
            minIndex = i;
            minDistance = distance;
        }
    }

    return minIndex;
}

ResultOfComparision Descriptors::compareDescriptors(const ListOfDesciptors& descriptors1,
                                                    const ListOfDesciptors& descriptors2){
    ResultOfComparision matches;

    for (unsigned i = 0; i < descriptors1.size(); i++) {
        auto first_min_index = getIndexOfNearest(descriptors1[i], descriptors2);
        auto second_min_index = getIndexOfNearest(descriptors2[first_min_index], descriptors1);
        if (second_min_index == i) {
            matches.emplace_back(i, first_min_index);
        }
    }

    return matches;
}


QImage Descriptors::getMergedMatrix(const Matrix& mat1,
                                    const Matrix& mat2,
                                    const Points& points1,
                                    const Points& points2,
                                    const ResultOfComparision& matches){
    Matrix mergedMat(std::max(mat1.getHeight(), mat2.getHeight()),
                          mat1.getWidth() + mat2.getWidth());
    for (int i = 0, ei = mat1.getHeight(); i < ei; i++) {
        for (int j = 0, ej = mat1.getWidth(); j < ej; j++) {
            mergedMat.setIntensity(i, j, mat1.getItensityAt(i, j));
        }
    }

    auto offset = mat1.getWidth();

    for (int i = 0, ei = mat2.getHeight(); i < ei; i++) {
        for (int j = 0, ej = mat2.getWidth(); j < ej; j++) {
            mergedMat.setIntensity(i, j + offset, mat2.getItensityAt(i, j));
        }
    }
    QImage result(Matrix::exportImage(mergedMat));

    QPainter painter(&result);
    painter.setPen(Qt::red);

    for (const auto& match : matches) {
        auto x1 = std::get<1>(points1[match.first]);
        auto y1 = std::get<0>(points1[match.first]);
        auto x2 = std::get<1>(points2[match.second]);
        auto y2 = std::get<0>(points2[match.second]);

        painter.drawEllipse(QPointF(x1, y1), 1, 1);
        painter.drawEllipse(QPointF(x2 + offset, y2), 1, 1);

        painter.drawLine(x1, y1, x2 + offset, y2);

    }

    return result;
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

double getMagnitude(const Desciptor& descriptor)
{
    double sum = 0;
    for(int i = 0; i < 128; i++)
    {
        sum += descriptor[i] * descriptor[i];
    }

    return sqrt(sum);
}

Desciptor Descriptors::Builder::normalize(const Desciptor& descriptor)
{
    Desciptor result = descriptor;
    double magnitude = getMagnitude(result);
    for(int i = 0; i < 128; i++) {
        result[i] /= magnitude;
        if (result[i] > 0.2) {
            result[i] = 0.2;
        }
    }

    magnitude = getMagnitude(result);
    for(int i = 0; i < 128; i++) {
        result[i] /= magnitude;
    }
    return result;
}

Descriptors::Builder& Descriptors::Builder::descriptors()
{
    printf("Building Descriptors\n");

    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    Matrix gradientValues = Sobel::Builder().sobelXY(sobelX, sobelY).build().getMatrix();
    Matrix gradientOrientations = Sobel::Builder().gradientOrientiations(sobelX, sobelY).build().getMatrix();

    const double binSize = M_PI*2 / this->numberOfBinsPerHistogram;

    for (const auto& point : this->filteredPoIs) {
        Desciptor descriptor(this->numberOfBins);
        //left-top corner point
        int x = std::get<0>(point) - this->gridCenter;
        int y = std::get<1>(point) - this->gridCenter;

        //building descriptor's histograms
        for(int i = 0; i < 16; i++){
            for(int j = 0; j < 16; j++){

                double gValue = gradientValues.getItensityAt(x + i, y + j);
                double gOrientation = gradientOrientations.getItensityAt(x + i, y + j);
                //indexing bin1
                int bin1Index = gOrientation / binSize;
                //check for end edge
                bin1Index %= this->numberOfBinsPerHistogram;
                double bin1Center = bin1Index * binSize + binSize / 2;

                //indexing bin2
                int bin2Index =  bin1Index + 1;
                if(gOrientation < bin1Center){
                    bin2Index  = bin1Index - 1;
                }
                //check for histogram's edges
                bin2Index = (bin2Index + this->numberOfBinsPerHistogram) % this->numberOfBinsPerHistogram;


                //get current histogram's index in grid
                int curHistogramIndexByX = i / this->histogramSize;
                int curHistogramIndexByY = j / this->histogramSize;

                int curHistogramIndex = curHistogramIndexByX  + curHistogramIndexByY * 4;

                //calculating distance to center
                double bin1Dist = abs(gOrientation - bin1Center);
                double bin2Dist = binSize - bin1Dist;

                descriptor[curHistogramIndex * numberOfBinsPerHistogram + bin1Index] +=
                        gValue * (1 - bin1Dist / binSize) ;
                descriptor[curHistogramIndex * numberOfBinsPerHistogram + bin2Index] +=
                        gValue * (1 - bin2Dist / binSize) ;
            }
        }

        //normalizing histograms
        descriptor = normalize(descriptor);

        this->listOfDesciptors.emplace_back(descriptor);
    }

    return *this;
}


Descriptors Descriptors::Builder::build() const
{
    return Descriptors(this->listOfDesciptors);
}
