#include <iostream>
//#include "src/model.h"
#include "KNHome.h"
#include "model.h"

using keras2cpp::Model;
using keras2cpp::Tensor;
using std::cout;
using std::endl;

void KNHome::testK2P() {
  // Initialize model.
  string fname = sparams["kerasModel"];

  auto kerasModel = Model::load(fname);

  // Create a 1D Tensor on length 10 for input data.
  Tensor in{2, 27};
  in.data_ = { {1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 4.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 3.0,  //0-26
               1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 4.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 3.0} }; //27-53

  cout << "model built by reading lib\n";
  in.print();
  // Run prediction.
  Tensor out = kerasModel(in);
  out.print();
  // cout << std::setprecision(8) << out(0, 0) << endl;
  in.data_[28] = 1.0;
  in.print();
  out = kerasModel(in);
  out.print();
  // cout << std::setprecision(8) << out(0, 0) << endl;

  return;
}
