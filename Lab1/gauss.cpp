#include "gauss.h"

Gauss::Gauss()
{

}

Pyramid Gauss::getPyramid() const
{
    return pyramid;
}

double Gauss::gauss(int x, int y, double sigma){
    double val = 2 * sigma * sigma;
    return 1.0 / (3.14 * val) * exp(-(x * x + y * y) / val);
}

Gauss::Builder::Builder(const Matrix& matrix)
    :matrix(matrix)
{}

Gauss::Builder& Gauss::Builder::init()
{
    printf("init\n");
    this->n = int(std::min(log2(matrix.getHeight() / minSize),
                      log2(matrix.getWidth() / minSize)));
    this->n++;
    printf("const: s = %d, sigmaA = %lf, sigma0 = %lf, n = %d\n", this->s, this->sigmaA, this->sigma0, n);
    return *this;
}

Gauss::Builder& Gauss::Builder::gaussPyramid(const Matrix& matrix){
    auto sigmaB = getSigmaB(this->sigma0, this->sigmaA);
    auto gauss = getGaussSeparable(sigmaB);

    this->pyramid.emplace_back (this->sigma0, this->sigma0, separableFilter(matrix, gauss));

    auto k = pow(2.0, 1.0/this->s);

    std::vector<Filter> filters;
    auto oldSigma = this->sigma0;
    for (int i = 0; i < this->s; i++) {
        auto newSigma = oldSigma * k;
        auto sigma = getSigmaB(newSigma, oldSigma);
        filters.emplace_back(getGaussSeparable(sigma));
        oldSigma = newSigma;
    }

    for (int oct = 0; oct < this->n; oct++) {
         for (int i = 1; i < this->s; i++) {
            Layers& layer = this->pyramid.back();
            double s0 = std::get<0>(layer);
            double s1 = std::get<1>(layer);

            Matrix& image = std::get<2>(layer);

            this->pyramid.emplace_back(s0*k, s1*k, separableFilter(image, filters[size_t(i)]));
         }

         if (oct != this->n - 1) {

             Layers& layer = this->pyramid.back();
             auto s = std::get<1>(layer);
             Matrix& image = std::get<2>(layer);

             this->pyramid.emplace_back(this->sigma0, s, Matrix::getDownscaled(image));
         }
    }
    return *this;
}

Filter Gauss::Builder::getGaussSeparable(const double sigmaB){
    int kernelSize = int(std::ceil(3 * sigmaB)) * 2 + 1;
    std::pair<Array,Array> separableFilter;
    separableFilter.first.resize(kernelSize);
    separableFilter.second.resize(kernelSize);

    int center = kernelSize / 2;

    for (int i = 0; i < kernelSize; i++) {
        double val =  sqrt(Gauss::gauss(i - center, i - center, sigmaB));
        separableFilter.first[i]= val;
        separableFilter.second[i]= val;
    }
    return separableFilter;

}

Matrix Gauss::Builder::separableFilter(const Matrix& matrix, const Filter& filter){
     Matrix tempMat(matrix.getHeight(),matrix.getWidth());
     int center = filter.first.size() / 2;
     for (int i = 0; i < matrix.getHeight(); i++) {
         for (int j = 0; j < matrix.getWidth(); j++) {
             double result_intensity = 0;
             for (std::size_t y = 0; y < filter.first.size(); y++) {
                 int jj = j + y - center;
                 jj = Matrix::getCol(jj,tempMat.getWidth());
                 auto intensity = matrix.getItensityAt(i, jj);
                 result_intensity += intensity * filter.first[y];
             }
             tempMat.setIntensity(i, j, result_intensity);
         }
     }

     Matrix result(tempMat.getHeight(), tempMat.getWidth());

     for (int i = 0; i < tempMat.getHeight(); i++) {
         for (int j = 0; j < tempMat.getWidth(); j++) {
             double result_intensity = 0;
             for (std::size_t x = 0; x < filter.second.size(); x++) {
                 int ii = i + x - center;
                 ii = Matrix::getRow(ii,tempMat.getHeight());
                 auto intensity = tempMat.getItensityAt(ii, j);
                 result_intensity += intensity * filter.second[x];
             }
             result.setIntensity(i, j, result_intensity);
         }
     }

     return result;
}

double Gauss::Builder::getSigmaB(double sigma0, double sigmaA) {
    return sqrt(sigma0 * sigma0 - sigmaA * sigmaA);
}

Gauss Gauss::Builder::build() const
{
    return Gauss(this->pyramid);
}
