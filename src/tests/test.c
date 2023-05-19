#include "leg2cheb.h"

int main(int argc, char *argv[]) {
  test_foreward_backward(1000, 100, 0, 1);
  test_foreward_backward(1000, 36, 1, 1);
  test_foreward_backward(10000, 100, 0, 1);
  test_foreward_backward(10000, 36, 1, 1);
  test_speed(100, 10, 100, 0, 18, 1);
  for (size_t i = 4; i < 12; i++)
    test_speed(36*pow(2, i+2), 36, 10, 0, 18, 1);
  test_direct(10, 0);
  test_direct(200, 0);
  test_direct(1000, 0);
  test_2_sizes(1000, 36, 1);
}
