#pragma once
#include <memory>
#include <QtGui>

class Array
{
public:

    Array(int size);

    double getItem(int index) const;

    void setItem(int index, double value);

    int getSize() const;
private:
    int size;
    std::vector<double> array;
};
