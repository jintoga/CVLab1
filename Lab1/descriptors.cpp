#include "descriptors.h"
#define PI 3.14
Descriptors::Descriptors()
{

}

ListOfDesciptors Descriptors::getDescriptors()
{
    return this->descriptors;
}
static double getDistance(const Desciptor& _first, const Desciptor& _second) {
    double sum = 0;
    for (size_t i = 0, ei = _first.size(); i < ei; i++) {
        sum += (_first[i] - _second[i])*(_first[i] - _second[i]);
    }
    return sqrt(sum);
}
static size_t getIndexOfNearest(const Desciptor& _descriptor, const ListOfDesciptors& _descriptors) {

    size_t min_index = 0;
    auto min_distance = getDistance(_descriptor, _descriptors[min_index]);

    for (size_t i = 1, ei = _descriptors.size(); i < ei; i++) {
        auto distance = getDistance(_descriptor, _descriptors[i]);
        if (distance < min_distance) {
            min_index = i;
            min_distance = distance;
        }
    }

    return min_index;
}

ResultOfComparision Descriptors::compareDescriptors(const ListOfDesciptors& descriptors1,
                                                    const ListOfDesciptors& descriptors2){
    ResultOfComparision matches;

    for (size_t i = 0, ei = descriptors1.size(); i < ei; i++) {
        auto first_min_index = getIndexOfNearest(descriptors1[i], descriptors2);
        auto second_min_index = getIndexOfNearest(descriptors2[first_min_index], descriptors1);
        if (second_min_index == i) {
            matches.emplace_back(i, first_min_index);
        }
    }

    return matches;
}
static void drawAroundPoint(QPainter& _painter, int _x, int _y) {
    _painter.drawPoint(_x - 1, _y);
    _painter.drawPoint(_x, _y - 1);
    _painter.drawPoint(_x, _y + 1);
    _painter.drawPoint(_x + 1, _y);
}
QImage Descriptors::getMergedMatrix(const Matrix& mat1,
                                    const Matrix& mat2,
                                    const Points& points1,
                                    const Points& points2,
                                    const ResultOfComparision& _matches){
    Matrix merged_image(std::max(mat1.getHeight(), mat2.getHeight()),
                          mat1.getWidth() + mat2.getWidth());
    for (int i = 0, ei = mat1.getHeight(); i < ei; i++) {
        for (int j = 0, ej = mat1.getWidth(); j < ej; j++) {
            merged_image.setIntensity(i, j, mat1.getItensityAt(i, j));
        }
    }

    auto offset = mat1.getWidth();

    for (int i = 0, ei = mat2.getHeight(); i < ei; i++) {
        for (int j = 0, ej = mat2.getWidth(); j < ej; j++) {
            merged_image.setIntensity(i, j + offset, mat2.getItensityAt(i, j));
        }
    }
    QImage result(Matrix::exportImage(merged_image));

    QPainter painter(&result);

    for (const auto& match : _matches) {
        auto x1 = std::get<1>(points1[match.first]);
        auto y1 = std::get<0>(points1[match.first]);
        auto x2 = std::get<1>(points2[match.second]);
        auto y2 = std::get<0>(points2[match.second]);

        int r = qrand() % 256, g = qrand() % 256, b = qrand() % 256;

        painter.setPen(QColor(r, g, b, 128));
        painter.drawLine(x1, y1, x2 + offset, y2);

        painter.setPen(QColor(r, g, b));
        drawAroundPoint(painter, x1, y1);
        drawAroundPoint(painter, x2 + offset, y2);
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
                if (bin1Index > numberOfBinsPerHistogram - 1) {
                    bin1Index = 0;
                    gOrientation -= M_PI * 2;
                }
                double bin1Center = bin1Index * binSize + binSize / 2;

                //indexing bin2
                int bin2Index = bin1Index + 1;
                if(gOrientation < bin1Center)
                    bin2Index  = bin1Index - 1;
                //check for histogram's edges
                bin2Index = (bin2Index + this->numberOfBinsPerHistogram) % this->numberOfBinsPerHistogram;


                //get current histogram's index in grid
                int curHistogramIndexByX = i / this->histogramSize;
                int curHistogramIndexByY = j / this->histogramSize;

                int curHistogramIndex = curHistogramIndexByX  + curHistogramIndexByY * 4;

                //calculating distance to center
                double bin1Dist = abs(gOrientation - bin1Center);

                double bin2Dist = binSize - bin1Dist;
                if (bin1Dist < 0 || bin2Dist < 0)
                    bin1Dist++;

                descriptor[curHistogramIndex * numberOfBinsPerHistogram + bin1Index] +=
                        gValue * (1 - bin1Dist / binSize) ;
                descriptor[curHistogramIndex * numberOfBinsPerHistogram + bin2Index] +=
                        gValue * (1 - bin2Dist / binSize) ;
            }
        }

        //normalizing histograms
        for(int i = 0; i < histogramSize*histogramSize; i++){
            double max = 1E-15;
            for(int j = 0; j < numberOfBinsPerHistogram; j++){
                if(max < descriptor[i*numberOfBinsPerHistogram + j]){
                    max = descriptor[i*numberOfBinsPerHistogram + j];
                }
            }
            for(int j = 0; j < numberOfBinsPerHistogram; j++){
                descriptor[i*numberOfBinsPerHistogram + j] = descriptor[i*numberOfBinsPerHistogram + j] / max;
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
