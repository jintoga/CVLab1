#include "array.h"

Array::Array(int size)
    :array(size)
{
}

double Array::getItem(int index) const {
    return array[index];
}

void Array::setItem(int index, double value){
    array[index] = value;
}

int Array::getSize() const {
    return array.capacity();
}
