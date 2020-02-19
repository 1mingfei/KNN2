#include "KNHome.h"


typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int me;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  auto t1 = Clock::now();
  KNHome* kn = new KNHome(argc, argv);
  if (kn) delete kn;
  auto t2 = Clock::now();
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0) {
    // ofstream ofs("log.txt", std::ofstream::app);
    cout << "Delta t2-t1: " << setprecision(8) \
         << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() \
         << " seconds" << endl;
  }

  MPI_Finalize();
  return (0);
}
