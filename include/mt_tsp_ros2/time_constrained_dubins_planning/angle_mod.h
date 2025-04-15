#include <math.h>

double fmodp(double x, double y) {
  double ret = fmod(x, y);
  if (ret < 0) {
    ret = y + ret;
  }
  return ret;
}

double angle_mod(double theta) {
  double twopi = 2*M_PI;
  double ret = fmodp(theta, twopi);
  while (ret > twopi) {
    ret -= twopi;
  }
  return ret;
}

double angdiff(double th1, double th2) {
  double th1_mod = angle_mod(th1);
  double th2_mod = angle_mod(th2);
  if (th2 - th1 > M_PI) {
    th2 -= 2*M_PI;
  } else if (th2 - th1 < -M_PI) {
    th2 += 2*M_PI;
  }
  return th2 - th1;
}
