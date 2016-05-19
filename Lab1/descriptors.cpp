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

Points Descriptors::getRPOIs()
{
    Points p(riDescriptors.size());
    for (size_t i = 0; i < p.size(); i++) {
        auto &d = riDescriptors[i];
        p[i] = make_tuple(get<1>(d), get<2>(d), get<3>(d));
    }
    return p;
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
        auto x1 = std::get<0>(points1[match.first]);
        auto y1 = std::get<1>(points1[match.first]);
        auto x2 = std::get<0>(points2[match.second]);
        auto y2 = std::get<1>(points2[match.second]);
        
        painter.drawEllipse(QPointF(x1, y1), 1, 1);
        painter.drawEllipse(QPointF(x2 + offset, y2), 1, 1);
        
        painter.drawLine(x1, y1, x2 + offset, y2);
        
    }
    
    return result;
}

QImage Descriptors::getMergedMatrix(const Matrix& mat1,
                                    const Matrix& mat2,
                                    const ListOfRIDescriptors& d1,
                                    const ListOfRIDescriptors& d2,
                                    const ResultOfComparision& _matches) {
    Points p1(d1.size());
    Points p2(d2.size());
    for (size_t i = 0; i < p1.size(); i++) {
        const auto& d = d1[i];
        p1[i] = make_tuple(get<1>(d), get<2>(d), get<3>(d));
    }
    for (size_t i = 0; i < p2.size(); i++) {
        const auto& d = d2[i];
        p2[i] = make_tuple(get<1>(d), get<2>(d), get<3>(d));
    }
    return Descriptors::getMergedMatrix(mat1, mat2, p1, p2, _matches);
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

Descriptors::Builder& Descriptors::Builder::descriptors()
{
    printf("Building Rotation Invariant Descriptors\n");
    const Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    const Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();
    Matrix gradientOrientations = Sobel::Builder().gradientOrientiations(sobelX, sobelY).build().getMatrix();
    
    const int DRAD = 8;
    double t = DRAD / 2;
    t *= t;
    Orientations orientations(binsOfWideHistogram);
    std::vector<int> peaks;
    
    
    for (const auto& point : this->filteredPoIs) {
        
        int x = std::get<0>(point);
        int y = std::get<1>(point);
        
        //init descriptor
        RorationInvariantDescriptor descriptor;
        std::get<1>(descriptor) = x;
        std::get<2>(descriptor) = y;
        Descriptor descriptorValue = std::get<0>(descriptor);
        descriptorValue.resize(numberOfBins);
        //clear orientations
        std::fill(begin(orientations), end(orientations), 0);
        
        //search for main orientation
        for (int i = -DRAD; i <= DRAD; ++i) {
            for (int j = -DRAD; j <= DRAD; ++j) {
                //skip if vector is bigger than descriptor's radius
                if (i * i + j * j > DRAD * DRAD)
                    continue;
                int qy = y + i;
                int qx = x + j;
                //check for edges of image
                if(qy < 0 || qy >= gradientOrientations.getHeight() 
                        || qx < 0 || qx >= gradientOrientations.getWidth()){
                    continue;
                }
                double dx = sobelX.getItensityAt(qy, qx);
                double dy = sobelY.getItensityAt(qy, qx);
                
                double len = sqrtl(dy * dy + dx * dx);
                //from Gauss formula
                len *= expl(-(i * i + j * j) / (2 * t));
                
                double fi = gradientOrientations.getItensityAt(qy, qx);
                double alph = fi * binsOfWideHistogram * 0.5 / M_PI;
                
                int bin1Index = int(alph) % binsOfWideHistogram;
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
        peaks.clear();
        int indexOfMax = 0;
        for (int i = 1; i < binsOfWideHistogram; i++) {
            if (orientations[i] > orientations[indexOfMax])
                indexOfMax = i;
        }
        peaks.push_back(indexOfMax);
        //2nd
        indexOfMax = (indexOfMax + 1) % binsOfWideHistogram;
        for (int i = 0; i < binsOfWideHistogram; i++) {
            if (i != peaks[0] && orientations[i] > orientations[indexOfMax])
                indexOfMax = i;
        }
        if (orientations[indexOfMax] >= orientations[peaks[0]] * 0.8)
            peaks.push_back(indexOfMax);
        
        
        for (unsigned i = 0; i < peaks.size(); i++) {
            //parabolic interpolation init
            int x2 = peaks[i];
            int x1 = (x2 + numberOfBins - 1) % numberOfBins;
            int x3 = (x2 + 1) % numberOfBins;
            double y1 = orientations[x1];
            double y2 = orientations[x2];
            double y3 = orientations[x3];
            
            //check if it is possible to make a parabol
            if (y2 <= y1 || y2 <= y3)
                continue;
            //interpolate
            double orientation = interpolation(x2, y1, y2, y3);
            get<3>(descriptor) = orientation;
            
            //building descriptor's histograms
            for (int i = -DRAD; i <= DRAD; i++) {
                for (int j = -DRAD; j <= DRAD; j++) {
                    //skip if vector is bigger than descriptor's radius
                    if (i * i + j * j > DRAD * DRAD)
                        continue;
                    int qy = y + i;
                    int qx = x + j;
                    //check for edges of image
                    if(qy < 0 || qy >= gradientOrientations.getHeight() 
                            || qx < 0 || qx >= gradientOrientations.getWidth()){
                        continue;
                    }
                    double dx = sobelX.getItensityAt(qy, qx);
                    double dy = sobelY.getItensityAt(qy, qx);
                    
                    double fi = gradientOrientations.getItensityAt(qy, qx) - orientation;
                    fi = fi < 0 ? fi + M_PI * 2 : fi; 
                    
                    double alph = fi * binsPerHistogram * 0.5 / M_PI;
                    
                    double len = sqrtf(dy * dy + dx * dx);
                    len *= expf(-(i * i + j * j) / (2 * t));
                    
                    int bin1Index = int(alph);
                    int bin2Index = bin1Index + 1;
                    if (bin2Index == binsPerHistogram)
                        bin2Index = 0;
                    
                    double weight = alph - bin1Index;
                    
                    //get normalized coordinates with vector rotation formula
                    qy = j * sinl(orientation) + i * cosl(orientation);
                    qx = j * cosl(orientation) - i * sinl(orientation);
                    
                    //get current histogram's index in grid
                    int curHistogramIndexByY = (qy + DRAD) * histogramSize / (DRAD * 2);
                    int curHistogramIndexByX = (qx + DRAD) * histogramSize / (DRAD * 2); 
                    
                    //check for rounding error
                    if (curHistogramIndexByY < 0)
                        curHistogramIndexByY = 0;
                    if (curHistogramIndexByY >= histogramSize)
                        curHistogramIndexByY = histogramSize - 1;
                    if (curHistogramIndexByX < 0)
                        curHistogramIndexByX = 0;
                    if (curHistogramIndexByX >= histogramSize)
                        curHistogramIndexByX = histogramSize - 1;   
                    
                    int curHistogramIndex = curHistogramIndexByX  + curHistogramIndexByY * histogramSize;
                    
                    descriptorValue[curHistogramIndex * binsPerHistogram + bin1Index] +=
                            len * (1 - weight);
                    descriptorValue[curHistogramIndex * binsPerHistogram + bin2Index] +=
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

