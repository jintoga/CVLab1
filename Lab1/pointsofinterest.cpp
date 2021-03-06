#include "pointsofinterest.h"
#include <cmath>

PointsOfInterest::PointsOfInterest()
{

}

Matrix PointsOfInterest::getMatrix(){
    return this->matrix;
}

Points PointsOfInterest::getPoIs(){
    return this->pois;
}

Points PointsOfInterest::getFilteredPoIs(){
    return this->filteredPoIs;
}

QImage PointsOfInterest::markPoints(const Matrix& mat, const Points& pois) {

    QImage result(Matrix::exportImage(mat));

    QPainter painter(&result);
    painter.setPen(Qt::red);

    for (const auto& point : pois) {
        auto y = std::get<0>(point);
        auto x = std::get<1>(point);
        painter.drawEllipse(QPointF(x,y), 1, 1);
    }

    return result;
}

PointsOfInterest::Builder::Builder()
{}

PointsOfInterest::Builder::Builder(const Matrix& matrix)
    :matrix(matrix)
{}

double PointsOfInterest::Builder::getContrast(const Matrix& matrix, const int w, const int x, const int y, const int dx, const int dy) {

    double sum = 0;
    auto center = w / 2;

    for (int i = -center; i <= center; i++) {
        for (int j = -center; j <= center; j++) {
            auto tmp =  matrix.getItensityAt(x + i, y + j) - matrix.getItensityAt(x + i + dx, y + j + dy);
            sum += tmp * tmp;
        }
    }

    return sum;
}

double PointsOfInterest::Builder::getDistance(int x1, int y1, int x2, int y2) {
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

PointsOfInterest::Builder& PointsOfInterest::Builder::init()
{
    printf("init\n");
    printf("const: w = %d, threshold = %lf, quantity = %d\n",
           this->w,
           this->threshold,
           this->quantity);
    return *this;
}

PointsOfInterest::Builder& PointsOfInterest::Builder::moravec()
{
    printf("Moravec\n");

    this->matrix = opMoravec(this->matrix, this->w);
    this->pois = findPoI(this->matrix, this->threshold);
    this->filteredPoIs = filterPoI(this->pois, this->quantity);


    return *this;
}

Matrix PointsOfInterest::Builder::opMoravec(const Matrix& matrix, int w) {

    Matrix result(matrix.getHeight(),matrix.getWidth());

    for (int i = 0; i < matrix.getHeight(); i++) {
        for (int j = 0; j < matrix.getWidth(); j++) {
            auto s = 99999999.0;        //value of operator
            //shift vector
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (dx != 0 && dy != 0) {
                        auto c = getContrast(matrix, w, i, j, dx, dy);     //contrast value of at point(i,j) in shift vector(dx, dy)
                        s = std::min(s, c);     //S(x,y) = min Contrast(x,y,d)
                    }
                }
            }
            result.setIntensity(i, j, s);
        }
    }

    return result;
}

PointsOfInterest::Builder& PointsOfInterest::Builder::harris()
{
    printf("Harris\n");

    this->matrix = opHarris(this->matrix, this->w);
    this->pois = findPoI(this->matrix, this->threshold);
    this->filteredPoIs = filterPoI(this->pois, this->quantity);


    return *this;
}

Matrix PointsOfInterest::Builder::opHarris(const Matrix& matrix, int w) {

    double k = 0.06;

    Matrix sobelX = Sobel::Builder(matrix).sobelX().build().getMatrix();
    Matrix sobelY = Sobel::Builder(matrix).sobelY().build().getMatrix();

    Matrix a(matrix.getHeight(),matrix.getWidth());
    Matrix b(matrix.getHeight(),matrix.getWidth());
    Matrix c(matrix.getHeight(),matrix.getWidth());

    auto center = w / 2;

    for (int i = 0; i < matrix.getHeight(); i++) {
        for (int j = 0; j < matrix.getWidth(); j++) {
            double sumA = 0, sumB = 0, sumC = 0;
            for (int x = -center; x <= center; x++) {
                for (int y = -center; y <= center; y++) {
                    auto ii = sobelX.getItensityAt(i + x, j + y);
                    auto jj = sobelY.getItensityAt(i + x, j + y);
                    sumA += ii * ii; //A = Ix*Ix
                    sumB += ii * jj; //B = Ix*Iy
                    sumC += jj * jj; //C = Iy*Iy
                }
            }
            a.setIntensity(i, j, sumA);
            b.setIntensity(i, j, sumB);
            c.setIntensity(i, j, sumC);
        }
    }

    Matrix result(matrix.getHeight(),matrix.getWidth());
    // f = det(h) - k*trace(H)2
    for (int i = 0; i < matrix.getHeight(); i++) {
        for (int j = 0; j < matrix.getWidth(); j++) {
            auto detH = a.getItensityAt(i, j) * c.getItensityAt(i, j) - b.getItensityAt(i, j)*b.getItensityAt(i, j);
            auto traceH = a.getItensityAt(i, j) + c.getItensityAt(i, j);
            auto harris = detH - k * traceH * traceH;
            result.setIntensity(i, j, harris);
        }
    }

    return result;
}

Points PointsOfInterest::Builder::findPoI(const Matrix& matrix, const double threshold) {

    Points points;

    for (int i = 0; i < matrix.getHeight(); i++) {
        for (int j = 0; j < matrix.getWidth(); j++) {
            if(matrix.getItensityAt(i, j) > threshold){
                bool isPoI = true;
                //run around neighbours
                for (int px = -1; px <= 1 && isPoI; px++) {
                    for (int py = -1; py <= 1 && isPoI; py++) {
                        if (px != 0 || py != 0) {
                            isPoI = matrix.getItensityAt(i, j) > matrix.getItensityAt(i + px, j + py);
                        }
                    }
                }
                if (isPoI) {
                    points.emplace_back(i, j, matrix.getItensityAt(i, j));
                }
            }

        }
    }

    return points;
}

Points PointsOfInterest::Builder::filterPoI(const Points& points, const unsigned quantity) {

    Points filteredPoI(points);

    int r = 0;

    while (filteredPoI.size() > quantity) {
        filteredPoI.erase(std::remove_if(filteredPoI.begin(), filteredPoI.end(),
                                         [&](const auto& _point) {
            for (const auto& point : filteredPoI) {
                auto distance = getDistance(std::get<0>(_point),
                                            std::get<1>(_point),
                                            std::get<0>(point),
                                            std::get<1>(point));
                if (distance < r && std::get<2>(_point) < std::get<2>(point)) {
                    return true;
                }
            }
            return false;
        }), filteredPoI.end());
        r++;
    }

    return filteredPoI;
}

PointsOfInterest PointsOfInterest::Builder::build() const
{
    return PointsOfInterest(this->matrix, this->pois, this->filteredPoIs);
}
