#include "descriptors.h"

Descriptors::Descriptors()
{

}

ListOfDescriptors Descriptors::getDescriptors()
{
    return this->descriptors;
}

static double getDistance(const Descriptor& desc1, const Descriptor& desc2) {
    double sum = 0;
    for (unsigned i = 0; i < desc1.size(); i++) {
        sum += (desc1[i] - desc2[i])*(desc1[i] - desc2[i]);
    }
    return sqrt(sum);
}

static unsigned getIndexOfNearest(const Descriptor& desc, const ListOfDescriptors& descriptors) {

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

ResultOfComparision Descriptors::compareDescriptors(const ListOfDescriptors& descriptors1,
                                                    const ListOfDescriptors& descriptors2){
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

double dist(const Descriptor& d)
{


    double sum = 0;
    for(int i = 0; i < 128; i++)
    {
        sum += (d[i] - d[i]) * (d[i] - d[i]);
    }

    return sqrt(sum);
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

Descriptor Descriptors::Builder::normalize(const Descriptor& descriptor)
{
    Descriptor result = descriptor;
    //normalizing histograms
    for(int i = 0; i < histogramSize*histogramSize; i++){
        double max = 1E-15;
        for(int j = 0; j < numberOfBinsPerHistogram; j++){
            if(max < result[i*numberOfBinsPerHistogram + j]){
                max = result[i*numberOfBinsPerHistogram + j];
            }
        }
        for(int j = 0; j < numberOfBinsPerHistogram; j++){
            result[i*numberOfBinsPerHistogram + j] = result[i*numberOfBinsPerHistogram + j] / max;
        }
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
        Descriptor descriptor(this->numberOfBins);
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
//                if (bin1Index > numberOfBinsPerHistogram - 1) {
//                    bin1Index = 0;
//                     gOrientation -= M_PI * 2;
//                }
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

        this->listOfDesciptors.emplace_back(normalize(descriptor    ));
    }

    return *this;
}

static int myProc(const int a, const int b)
{
    return (a + b) % b;
}

static double interpolate(const double x0, const double x1, const double x2,
                   const double y0, const double y1, const double y2)
{
    double b = (y2 - y0) / (x2 - x0) / (x2 - x1) - (y1 - y0) / (x1 - x0) / (x2 - x1);
    double a = (y1 - y0) / (x1 - x0) - b * (x1 + x0);

    return -a / 2 / b;
}

std::pair<double, double> rotate(double x0, double y0, double x1, double y1, double angle)
{
    if(x0 == x1 && y0 == y1)
        return std::make_pair(x0, y0);

    //находим угол
    double deltaX = x1 - x0;
    double deltaY = y0 - y1;
    double dist = hypot(deltaX, deltaY);

    double angleQQ = acos(deltaX / dist);
    if(deltaY < 0)
        angleQQ = (2 * M_PI - angleQQ) / M_PI * 180;
    else
        angleQQ = angleQQ / M_PI * 180;

    double newAngle = angleQQ - angle;
    if(newAngle < 0)
        newAngle += 360;

    double newDeltaX = cos(newAngle / 180 * M_PI) * dist;
    double newDeltaY = sin(newAngle / 180 * M_PI) * dist;

    double newX = x0 + newDeltaX;
    double newY = y0 - newDeltaY;


//    double newX = (x0 - 8) * cos(angle) - (y0 - 8) * sin(angle);
//    double newY = (x1 - 8) * sin(angle) + (y1 - 8) * cos(angle);

    return std::make_pair(newX,newY);
}

std::pair<int, int> findMaxPair(const Orientations& orientations)
{
    int res1 = 0, res2 = 1;
    double max1 = orientations[0];
    double max2 = orientations[1];
    if(max1 < max2)
    {
        std::swap(max1, max2);
        std::swap(res1, res2);
    }


    for(unsigned i = 2; i < orientations.size(); i++)
    {
        double val = orientations[i];
        if(val > max1)     //нашли меньше обоих
        {
            std::swap(max1, max2);
            std::swap(res1, res2);
            max1 = val;
            res1 = i;
        }
        else
            if(val > max2)     //нашли меньше второго
            {
                max2 = val;
                res2 = i;
            }
    }
    return std::make_pair(res1, res2);
}

Orientations Descriptors::Builder::findOrientationBins(int x, int y, double orientationBinSize)
{
    Orientations orientations(this->numberOfOrientationBins);
    for(int i = 0; i < 16; i++){
        for(int j = 0; j < 16; j++){
            double gValue = gradientValues.getItensityAt(x + i, y + j);
            double gOrientation = gradientOrientations.getItensityAt(x + i, y + j);

            //indexing bin1
            int bin1Index = gOrientation / orientationBinSize;
            //check for end edge
            bin1Index %= this->numberOfOrientationBins;
            if (bin1Index > numberOfOrientationBins - 1) {
                bin1Index = 0;
                gOrientation -= M_PI * 2;
            }
            double bin1Center = bin1Index * orientationBinSize + orientationBinSize / 2;

            //indexing bin2
            int bin2Index =  bin1Index + 1;
            if(gOrientation < bin1Center){
                bin2Index  = bin1Index - 1;
            }
            //check for histogram's edges
            bin2Index = (bin2Index + this->numberOfOrientationBins) % this->numberOfOrientationBins;

            //calculating distance to center
            double bin1Dist = abs(gOrientation - bin1Center);
            double bin2Dist = orientationBinSize - bin1Dist;


            orientations[bin1Index] += gValue * (1 - bin1Dist / orientationBinSize) ;
            orientations[bin2Index] += gValue * (1 - bin2Dist / orientationBinSize) ;

        }
    }

    return orientations;
}

Descriptor Descriptors::Builder::getFinalBins(Point point, double mainOrt, int x, int y, double binSize)
{
    Descriptor descriptor(this->numberOfBins);

    //max distance
    int maxDist = this->gridCenter * sqrt(2) + 1;
    int pX = std::get<0>(point);
    int pY = std::get<1>(point);
    for(int i = pX - maxDist; i < pX; i++){
        for(int j = pY - maxDist; j < pY; j++){
            auto temp = rotate(pX, pY, i, j, mainOrt);
            double newX = temp.first;
            double newY = temp.second;
            //after rotation check if new coordinates are inside of current area
            if(newX < x || newX > x + 16 || newY < y || newY > y + 16){
                continue;
            }

            newX -= x;
            newY -= y;

            //get current histogram's index in grid
            int curHistogramIndexByX = newX / this->histogramSize;
            int curHistogramIndexByY = newY / this->histogramSize;

            int curHistogramIndex = curHistogramIndexByX  + curHistogramIndexByY * 4;

            if(curHistogramIndex >31 || curHistogramIndex <0)
                curHistogramIndex = 1;

            double gValue = gradientValues.getItensityAt(i, j);
            double gOrientation = gradientOrientations.getItensityAt(i, j) - mainOrt;
            if(gOrientation < 0)
                gOrientation += 2 * M_PI;

            //indexing bin1
            int bin1Index = gOrientation / binSize;
            //check for end edge
            bin1Index %= this->numberOfBinsPerHistogram;
            if (bin1Index > numberOfBinsPerHistogram - 1) {
                bin1Index = 0;
                gOrientation -= M_PI * 2;
            }
            double bin1Center = bin1Index * binSize + binSize / 2;

            //indexing bin2
            int bin2Index =  bin1Index + 1;
            if(gOrientation < bin1Center){
                bin2Index  = bin1Index - 1;
            }
            //check for histogram's edges
            bin2Index = (bin2Index + this->numberOfBinsPerHistogram) % this->numberOfBinsPerHistogram;

            //calculating distance to center
            double bin1Dist = abs(gOrientation - bin1Center);
            double bin2Dist = binSize - bin1Dist;

            descriptor[curHistogramIndex * numberOfBinsPerHistogram + bin1Index] +=
                    gValue * (1 - bin1Dist / binSize) ;
            descriptor[curHistogramIndex * numberOfBinsPerHistogram + bin2Index] +=
                    gValue * (1 - bin2Dist / binSize) ;
        }
    }


    return descriptor;
}

Descriptors::Builder& Descriptors::Builder::invariantRotationDescriptors()
{
    printf("Building Invariant Rotation Descriptors\n");

    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    this->gradientValues = Sobel::Builder().sobelXY(sobelX, sobelY).build().getMatrix();
    this->gradientOrientations = Sobel::Builder().gradientOrientiations(sobelX, sobelY).build().getMatrix();

    const double binSize = M_PI*2 / this->numberOfBinsPerHistogram;
    const double orientationBinSize = M_PI*2 / this->numberOfOrientationBins;
    for (const auto& point : this->filteredPoIs) {

        //left-top corner point
        int x = std::get<0>(point) - this->gridCenter;
        int y = std::get<1>(point) - this->gridCenter;

        Orientations orientations = findOrientationBins(x, y, orientationBinSize);

        //find main bins
       auto mainBins = findMaxPair(orientations);
       //indexing bin1
       int bin1Index = mainBins.first;
       //find main orientation
       double orientation0 = bin1Index * orientationBinSize - orientationBinSize / 2;
       double orientation1 = bin1Index * orientationBinSize + orientationBinSize / 2;
       double orientation2 = bin1Index * orientationBinSize + orientationBinSize * 3 / 2;

       double mainOrt = interpolate(orientation0, orientation1, orientation2,
                                    orientations[myProc(bin1Index - 1,numberOfBinsPerHistogram)],
                                    orientations[bin1Index],
                                    orientations[myProc(bin1Index + 1,numberOfBinsPerHistogram)]);

       //final distribution
       Descriptor descriptor(this->numberOfBins);
       descriptor = getFinalBins(point, mainOrt, x, y, binSize);
       this->listOfDesciptors.emplace_back(normalize(descriptor));

       //check if there is the 2nd peak
       int bin2Index = mainBins.second;
       if(orientations[bin2Index] > 0.8 * orientations[bin1Index]){
           double orientation0 = bin2Index * orientationBinSize - orientationBinSize / 2;
           double orientation1 = bin2Index * orientationBinSize + orientationBinSize / 2;
           double orientation2 = bin2Index * orientationBinSize + orientationBinSize * 3 / 2;
           double mainOrt = interpolate(orientation0, orientation1, orientation2,
                                        orientations[myProc(bin2Index - 1,numberOfBinsPerHistogram)],
                                        orientations[bin2Index],
                                        orientations[myProc(bin2Index + 1,numberOfBinsPerHistogram)]);
           //final distribution
           Descriptor descriptor(this->numberOfBins);
           descriptor = getFinalBins(point, mainOrt, x, y, binSize);
           this->listOfDesciptors.emplace_back(normalize(descriptor));
       }
    }

    return *this;
}


Descriptors Descriptors::Builder::build() const
{
    return Descriptors(this->listOfDesciptors);
}
