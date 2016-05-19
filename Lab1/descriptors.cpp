#include "descriptors.h"

Descriptors::Descriptors()
{

}

ListOfDescriptors Descriptors::getDescriptors(boolean isRotationInvariant)
{
    if(isRotationInvariant){
        ListOfDescriptors listOfDescriptors;
        listOfDescriptors.reserve(this->riDescriptors.size());
        for(RorationInvariantDescriptor riDesc : riDescriptors){
            listOfDescriptors.emplace_back(std::get<0>(riDesc));
        }
        return listOfDescriptors;
    }
    return this->descriptors;
}

ListOfRIDescriptors Descriptors::getRIDescriptors()
{
    return this->riDescriptors;
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

ResultOfComparision Descriptors::compareDescriptors(const ListOfRIDescriptors& descriptors1,
                                                    const ListOfRIDescriptors& descriptors2){
    ResultOfComparision matches;
    matches.reserve(min(descriptors1.size(), descriptors2.size()));
    vector<size_t> lp(descriptors1.size());
    vector<size_t> rp(descriptors2.size());

    vector<double> left(descriptors1.size());
    vector<double> right(descriptors2.size());
    fill(begin(left), end(left), numeric_limits<double>::max());
    fill(begin(right), end(right), numeric_limits<double>::max());

    double thres = 7e-1;
    for (size_t i = 0; i < descriptors1.size(); i++) {
        for (size_t j = 0; j < descriptors2.size(); j++) {
            double cur = 0;
            for (size_t k = 0; k < 128; k++) {
                double dim = get<0>(descriptors1[i])[k] - get<0>(descriptors2[j])[k];
                cur += dim * dim;
            }
            if (left[i] > cur) {
                left[i] = cur;
                lp[i] = j;
            }
            if (right[j] > cur) {
                right[j] = cur;
                rp[j] = i;
            }
        }
    }

    for (size_t i = 0; i < descriptors1.size(); i++) {
        int j = lp[i];
        if (rp[j] == i && left[i] <= thres * thres) {
            matches.push_back(make_pair(i, j));
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

double getMagnitude(const Descriptor& descriptor, int numberOfBins)
{
    double sum = 0;
    for(int i = 0; i < numberOfBins; i++)
    {
        sum += descriptor[i] * descriptor[i];
    }

    return sqrt(sum);
}

Descriptor Descriptors::Builder::normalize(const Descriptor& descriptor)
{
    Descriptor result = descriptor;
    double magnitude = getMagnitude(result, numberOfBins);
    for(int i = 0; i < numberOfBins; i++) {
        result[i] /= magnitude;
        if (result[i] > 0.2) {
            result[i] = 0.2;
        }
    }

    magnitude = getMagnitude(result, numberOfBins);
    for(int i = 0; i < numberOfBins; i++) {
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

    const double binSize = M_PI*2 / this->binsPerHistogram;

    for (const auto& point : this->filteredPoIs) {
        Descriptor descriptor(this->numberOfBins);
        //left-top corner point
        int x = std::get<0>(point) - this->gridCenter;
        int y = std::get<1>(point) - this->gridCenter;

        //building descriptor's histograms
        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16; j++) {

                double gValue = gradientValues.getItensityAt(x + i, y + j);
                double gOrientation = gradientOrientations.getItensityAt(x + i, y + j);
                //indexing bin1
                int bin1Index = gOrientation / binSize;
                //check for end edge
                bin1Index %= this->binsPerHistogram;
                double bin1Center = bin1Index * binSize + binSize / 2;

                //indexing bin2
                int bin2Index =  bin1Index + 1;
                if(gOrientation < bin1Center){
                    bin2Index  = bin1Index - 1;
                }
                //check for histogram's edges
                bin2Index = (bin2Index + this->binsPerHistogram) % this->binsPerHistogram;


                //get current histogram's index in grid
                int curHistogramIndexByX = i / this->histogramSize;
                int curHistogramIndexByY = j / this->histogramSize;

                int curHistogramIndex = curHistogramIndexByX  + curHistogramIndexByY * 4;

                //calculating distance to center
                double bin1Dist = abs(gOrientation - bin1Center);
                double bin2Dist = binSize - bin1Dist;

                descriptor[curHistogramIndex * binsPerHistogram + bin1Index] +=
                        gValue * (1 - bin1Dist / binSize) ;
                descriptor[curHistogramIndex * binsPerHistogram + bin2Index] +=
                        gValue * (1 - bin2Dist / binSize) ;
            }
        }

        //normalizing histograms
        descriptor = normalize(descriptor);

        this->listOfDescriptors.emplace_back(descriptor);
    }

    return *this;
}

double Descriptors::Builder::interpolation(const double x2, const double y1, const double y2, const double y3)
{
        double ax[] = {double(x2 - 1), double(x2), double(x2 + 1)};
        double ay[] = {y1, y2, y3};

        double q1[3] = {y1, 0, 0}; // current coeffs
        double q2[3]; // additional coeffs

        //parabolic interpolation
        for (int i = 1; i < 3; i++) {
            double co = 1;
            for (int j = 0; j < i; j++) {
                co *= (ax[i] - ax[j]);
            }
            double cur = ax[i] * (ax[i] * q1[2] + q1[1]) + q1[0];
            co = (ay[i] - cur) / co;
            fill(begin(q2), end(q2), 0);
            q2[0] = 1;
            for (int j = 0; j < i; j++) {
                double z = 0;
                for (int c3 = 0; c3 <= j; c3++) {
                    double nx = q2[c3];
                    q2[c3] = z - (nx * ax[j]);
                    z = nx;
                }
                q2[j + 1] = z;
            }
            for (int j = 0; j < 3; j++) {
                q1[j] += q2[j] * co;
            }
        }

        // resulting orientation
        double res = -q1[1] / (2 * q1[2]);
        if (res < 0)
            res += binsOfWideHistogram;
        return M_PI * 2 * (res) / binsOfWideHistogram;
}

Descriptors::Builder& Descriptors::Builder::rotationInvariantDescriptors()
{
    printf("Building Rotation Invariant Descriptors\n");
    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();
    Matrix gradientOrientations = Sobel::Builder().gradientOrientiations(sobelX, sobelY).build().getMatrix();

    const int DRAD = 4;
    double sqs = DRAD / 2;
    sqs *= sqs;
    Orientations orientations(binsOfWideHistogram);
    std::vector<int> dirs;
    for (const auto& point : this->filteredPoIs) {

        int x = std::get<0>(point) + DRAD;
        int y = std::get<1>(point) + DRAD;

        //init descriptor
        RorationInvariantDescriptor descriptor;
        std::get<1>(descriptor) = x - DRAD;
        std::get<2>(descriptor) = y - DRAD;
        Descriptor descriptorValue = std::get<0>(descriptor);
        descriptorValue.resize(numberOfBins);
        //clear orientations
        std::fill(begin(orientations), end(orientations), 0);

        //search for main orientation
        for (int cy = -DRAD; cy < DRAD; ++cy) {
            for (int cx = -DRAD; cx < DRAD; ++cx) {
                if (cy * cy + cx * cx > DRAD * DRAD)
                    continue;
                int qy = x + cy;
                int qx = y + cx;
                double dx = sobelX.getItensityAt(qx, qy);
                double dy = sobelY.getItensityAt(qx, qy);

                double len = sqrtf(dy * dy + dx * dx);
                len *= expf(-(cy * cy + cx * cx) / (2.f * sqs));

                double fi = gradientOrientations.getItensityAt(qx, qy);
                double alph = fi * binsOfWideHistogram * 0.5f / M_PI;

                int bin1Index = int(alph);
                int bin2Index = bin1Index + 1;
                if (bin2Index == binsOfWideHistogram)
                    bin2Index = 0;

                double weight = alph - bin1Index;
                orientations[bin1Index] += len * (1 - weight);
                orientations[bin2Index] += len * weight;
            }
        }

        //find max pair
        //1st
        dirs.clear();
        int indexOfMax = 0;
        for (int i = 1; i < binsOfWideHistogram; i++) {
            if (orientations[i] > orientations[indexOfMax])
                indexOfMax = i;
        }
        dirs.push_back(indexOfMax);
        //2nd
        indexOfMax = (indexOfMax + 1) % binsOfWideHistogram;
        for (int i = 0; i < binsOfWideHistogram; i++) {
            if (i != dirs[0] && orientations[i] > orientations[indexOfMax])
                indexOfMax = i;
        }
        if (orientations[indexOfMax] >= orientations[dirs[0]] * 0.8)
            dirs.push_back(indexOfMax);


        for (unsigned i = 0; i < dirs.size(); i++) {
            // interpolation init
            int x2 = dirs[i];
            int x1 = (x2 + numberOfBins - 1) % numberOfBins;
            int x3 = (x2 + 1) % numberOfBins;
            double y1 = orientations[x1];
            double y2 = orientations[x2];
            double y3 = orientations[x3];

            if (y2 < y1 || y2 < y3)
                continue;
            //interpolate
            double orientation = interpolation(x2, y1, y2, y3);
            get<3>(descriptor) = orientation;

            double rsin = sinf(orientation);
            double rcos = cosf(orientation);

            for (int cy = -DRAD; cy <= DRAD; cy++) {
                for (int cx = -DRAD; cx <= DRAD; cx++) {
                    if (cy * cy + cx * cx > DRAD * DRAD)
                        continue;

                    double dx = sobelX.getItensityAt(y + cx, x + cy);
                    double dy = sobelY.getItensityAt(y + cx, x + cy);

                    double fi = gradientOrientations.getItensityAt(y + cx, x + cy) - orientation;
                    if (fi < 0)
                        fi += M_PI * 2.;
                    double alph = fi * binsPerHistogram * .5 / M_PI;

                    double len = sqrtf(dy * dy + dx * dx);
                    len *= expf(-(cy * cy + cx * cx) / (2.f * sqs));

                    int bin1Index = int(alph);
                    int bin2Index = bin1Index + 1;
                    if (bin2Index == binsPerHistogram)
                        bin2Index = 0;

                    double weight = alph - bin1Index;

                    // getting normalized coorditates
                    double qy = -((-cy) * rcos - (cx) * rsin);
                    double qx = (cx) * rcos + (-cy) * rsin;
                    // box`s selection
                    int ybox = int((qy - 0.5 + DRAD) / histogramSize);
                    int xbox = int((qx - 0.5 + DRAD) / histogramSize);
                    if (ybox < 0)
                        ybox = 0;
                    if (ybox >= histogramSize)
                        ybox = histogramSize - 1;
                    if (xbox < 0)
                        xbox = 0;
                    if (xbox >= histogramSize)
                        xbox = histogramSize - 1;

                    descriptorValue[(ybox * histogramSize + xbox) * binsPerHistogram + bin1Index] +=
                            len * (1 - weight);
                    descriptorValue[(ybox * histogramSize + xbox) * binsPerHistogram + bin2Index] +=
                            len * weight;
                }
            }
            std::get<0>(descriptor) = normalize(descriptorValue);
            this->listOfRIDescriptors.emplace_back(descriptor);

        }

    }

    return *this;
}

Descriptors Descriptors::Builder::build() const
{
    return Descriptors(this->listOfRIDescriptors);
}
