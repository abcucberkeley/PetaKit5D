/*
 * replace_nan_with_value.cpp - replace nan in matrix
 *
 * Implemented by Xiongtao Ruan.
 * 
 * This is a MEX-file for MATLAB.
 */ 

#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"


using namespace matlab::data;
using matlab::mex::ArgumentList;


using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        TypedArray<double> largeArray = std::move(inputs[0]);
        TypedArray<double> givenVal = inputs[1];
        for (uint64_t i = 0; i < largeArray.getNumberOfElements(); i++) {
            if (std::isnan(*(largeArray[i]))) {
                largeArray[i] = givenVal[0];
            }
        }
        outputs[0] = largeArray;
    }
};

