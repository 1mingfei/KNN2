#include <iostream>
#include "model.h"

using keras2cpp::Model;
using keras2cpp::Tensor;
using std::cout;
using std::endl;

int main() {
    // Initialize model.
    auto model = Model::load("keras.model.1");

    // Create a 1D Tensor on length 10 for input data.
    Tensor in{2, 27};
    in.data_ = { {1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 4.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 3.0,  //0-26
                 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 4.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 3.0} }; //27-53
    in.print();
    // Run prediction.
    Tensor out = model(in);
    out.print();
    // cout << std::setprecision(8) << out(0, 0) << endl;
    in.data_[28] = 1.0;
    in.print();
    out = model(in);
    out.print();
    // cout << std::setprecision(8) << out(0, 0) << endl;

    return 0;
}
