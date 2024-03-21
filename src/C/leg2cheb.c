#include "leg2cheb.h"

size_t lmin(size_t a, size_t b) { return (a > b) ? b : a; }

inline static double chebval(const double x, const double *c, const size_t N) {
  const double x2 = 2 * x;
  double c0 = c[N - 2];
  double c1 = c[N - 1];
  for (size_t i = 3; i < N + 1; i++) {
    double tmp = c0;
    c0 = c[N - i] - c1;
    c1 = tmp + c1 * x2;
  }
  return c0 + c1 * x;
}

void __Lambda(const double *z, double *w, const size_t N) {
  for (size_t i = 0; i < N; i++)
    w[i] = _Lambda(z[i]);
}

double _Lambda0(const double z) { return exp(lgamma(z + 0.5) - lgamma(z + 1)); }

inline static double _LambdaE10(const double z) {
  const double a0[6] = {9.9992197073403921e-01, -7.7997420312031631e-05,
                        3.1806372211417961e-08, -3.9175710357440109e-11,
                        1.0041660721179864e-13, -4.3833711744079551e-16};
  const double z0 = z + 0.25;
  return chebval(-1 + 200 / pow(z0, 2), a0, 6) / sqrt(z0);
}

inline static double _LambdaE24(const double z) {
  const double a0[4] = {9.9998643952730293e-01, -1.3559507925471981e-05,
                        9.6456296275015755e-10, -2.0852410532893096e-13};
  const double z0 = z + 0.25;
  return chebval(-1 + 1152 / pow(z0, 2), a0, 4) / sqrt(z0);
}

inline static double _LambdaE64(const double z) {
  const double a0[3] = {9.9999809270865958e-01, -1.9072722439854597e-06,
                        1.9095897782497945e-11};
  const double z0 = z + 0.25;
  return chebval(-1 + 8192 / pow(z0, 2), a0, 3) / sqrt(z0);
}

inline static double _LambdaE600(const double z) {
  const double a0[2] = {9.9999997829861853e-01, -2.1701378998945873e-08};
  const double z0 = z + 0.25;
  // return chebval(-1 + 720000 / pow(z0, 2), a0, 2) / sqrt(z0);
  return (a0[0] + a0[1] * (-1 + 720000 / pow(z0, 2))) / sqrt(z0);
}

double _LambdaE(const double z) {
  if (z > 600) {
    return _LambdaE600(z);
  } else if (z > 64) {
    return _LambdaE64(z);
  } else if (z > 24) {
    return _LambdaE24(z);
  }
  return _LambdaE10(z);
}

double _Lambda(const double z) {
  double z2, zz, y, s;
  double zp = z + 0.25;
  z2 = zp * zp;
  s = sqrt(zp);
  y = 1 - 0.015625 / z2;
  if (z > 800)
    return y / s;
  zz = z2 * z2;
  y += 0.0025634765625 / zz;
  if (z > 90)
    return y / s;
  zz *= z2;
  y -= 0.0012798309326171875 / zz;
  if (z > 32)
    return y / s;
  zz *= z2;
  y += 0.0013435110449790955 / zz;
  if (z > 17)
    return y / s;
  zz *= z2;
  y -= 0.0024328966392204165 / zz;
  if (z > 12)
    return y / s;
  zz *= z2;
  y += 0.006754237533641572 / zz;
  if (z > 10)
    return y / s;
  return _Lambda0(z);
}

const double LamInt[64] = {
    1.7724538509055161e+00, 8.8622692545275805e-01, 6.6467019408956851e-01,
    5.5389182840797380e-01, 4.8465534985697706e-01, 4.3618981487127934e-01,
    3.9984066363200604e-01, 3.7128061622971992e-01, 3.4807557771536241e-01,
    3.2873804562006448e-01, 3.1230114333906128e-01, 2.9810563682364938e-01,
    2.8568456862266400e-01, 2.7469670059871537e-01, 2.6488610414876129e-01,
    2.5605656734380255e-01, 2.4805479961430874e-01, 2.4075907021388790e-01,
    2.3407131826350211e-01, 2.2791154673025205e-01, 2.2221375806199575e-01,
    2.1692295429861491e-01, 2.1199288715546458e-01, 2.0738434613034576e-01,
    2.0306383891929691e-01, 1.9900256214091097e-01, 1.9517558979204730e-01,
    1.9156122701812048e-01, 1.8814049082136833e-01, 1.8489668925548267e-01,
    1.8181507776789130e-01, 1.7888257651357048e-01, 1.7608753625554593e-01,
    1.7341954328197706e-01, 1.7086925588077151e-01, 1.6842826651104620e-01,
    1.6608898503172612e-01, 1.6384453928805415e-01, 1.6168869008689554e-01,
    1.5961575816270457e-01, 1.5762056118567075e-01, 1.5569835921999184e-01,
    1.5384480732451575e-01, 1.5205591421609116e-01, 1.5032800609999922e-01,
    1.4865769492111033e-01, 1.4704185041109827e-01, 1.4547757540672487e-01,
    1.4396218399623817e-01, 1.4249318211872553e-01, 1.4106825029753828e-01,
    1.3968522823579768e-01, 1.3834210104122271e-01, 1.3703698688045646e-01,
    1.3576812589082260e-01, 1.3453387020090604e-01, 1.3333267493125509e-01,
    1.3216309006343707e-01, 1.3102375308013156e-01, 1.2991338229131691e-01,
    1.2883077077222260e-01, 1.2777478084786012e-01, 1.2674433906682897e-01,
    1.2573843161391765e-01};

// Lambda with integer argument. For initializing direct solver
double _LambdaI(const size_t z) {
  if (z > 63)
    return _LambdaE64((double)z);
  return LamInt[z];
}

// Row major DCT II matrix of shape 18 x 18
const double DCTII[324] = {
    5.5555555555555552e-02,  5.5555555555555552e-02,  5.5555555555555552e-02,
    5.5555555555555552e-02,  5.5555555555555552e-02,  5.5555555555555552e-02,
    5.5555555555555552e-02,  5.5555555555555552e-02,  5.5555555555555552e-02,
    5.5555555555555552e-02,  5.5555555555555552e-02,  5.5555555555555552e-02,
    5.5555555555555552e-02,  5.5555555555555552e-02,  5.5555555555555552e-02,
    5.5555555555555552e-02,  5.5555555555555552e-02,  5.5555555555555552e-02,
    1.1068829978797173e-01,  1.0732509180989648e-01,  1.0070086522629444e-01,
    9.1016893809887978e-02,  7.8567420131838622e-02,  6.3730715150116246e-02,
    4.6957584637855494e-02,  2.8757671678057886e-02,  9.6839714164064592e-03,
    -9.6839714164064453e-03, -2.8757671678057872e-02, -4.6957584637855480e-02,
    -6.3730715150116246e-02, -7.8567420131838608e-02, -9.1016893809887991e-02,
    -1.0070086522629444e-01, -1.0732509180989648e-01, -1.1068829978797173e-01,
    1.0942308366802311e-01,  9.6225044864937631e-02,  7.1420845520726597e-02,
    3.8002238147296537e-02,  6.8035933285964064e-18,  -3.8002238147296502e-02,
    -7.1420845520726597e-02, -9.6225044864937603e-02, -1.0942308366802311e-01,
    -1.0942308366802313e-01, -9.6225044864937617e-02, -7.1420845520726611e-02,
    -3.8002238147296502e-02, -2.0410779985789219e-17, 3.8002238147296558e-02,
    7.1420845520726584e-02,  9.6225044864937645e-02,  1.0942308366802311e-01,
    1.0732509180989648e-01,  7.8567420131838622e-02,  2.8757671678057886e-02,
    -2.8757671678057872e-02, -7.8567420131838608e-02, -1.0732509180989648e-01,
    -1.0732509180989649e-01, -7.8567420131838636e-02, -2.8757671678057848e-02,
    2.8757671678057810e-02,  7.8567420131838595e-02,  1.0732509180989648e-01,
    1.0732509180989649e-01,  7.8567420131838636e-02,  2.8757671678057959e-02,
    -2.8757671678057799e-02, -7.8567420131838650e-02, -1.0732509180989648e-01,
    1.0441029119843427e-01,  5.5555555555555566e-02,  -1.9294241962992256e-02,
    -8.5116049235441985e-02, -1.1111111111111110e-01, -8.5116049235442040e-02,
    -1.9294241962992259e-02, 5.5555555555555483e-02,  1.0441029119843427e-01,
    1.0441029119843430e-01,  5.5555555555555525e-02,  -1.9294241962992207e-02,
    -8.5116049235442040e-02, -1.1111111111111110e-01, -8.5116049235441971e-02,
    -1.9294241962992287e-02, 5.5555555555555629e-02,  1.0441029119843426e-01,
    1.0070086522629444e-01,  2.8757671678057886e-02,  -6.3730715150116246e-02,
    -1.1068829978797173e-01, -7.8567420131838636e-02, 9.6839714164064315e-03,
    9.1016893809887950e-02,  1.0732509180989649e-01,  4.6957584637855473e-02,
    -4.6957584637855411e-02, -1.0732509180989648e-01, -9.1016893809887922e-02,
    -9.6839714164064991e-03, 7.8567420131838511e-02,  1.1068829978797173e-01,
    6.3730715150116218e-02,  -2.8757671678057581e-02, -1.0070086522629434e-01,
    9.6225044864937631e-02,  6.8035933285964064e-18,  -9.6225044864937603e-02,
    -9.6225044864937617e-02, -2.0410779985789219e-17, 9.6225044864937645e-02,
    9.6225044864937673e-02,  3.4017966642982035e-17,  -9.6225044864937645e-02,
    -9.6225044864937687e-02, -4.7625153300174848e-17, 9.6225044864937631e-02,
    9.6225044864937687e-02,  6.1232339957367660e-17,  -9.6225044864937520e-02,
    -9.6225044864937700e-02, 1.2253345554102288e-16,  9.6225044864937617e-02,
    9.1016893809887978e-02,  -2.8757671678057872e-02, -1.1068829978797173e-01,
    -4.6957584637855550e-02, 7.8567420131838595e-02,  1.0070086522629444e-01,
    -9.6839714164064176e-03, -1.0732509180989648e-01, -6.3730715150116357e-02,
    6.3730715150116135e-02,  1.0732509180989644e-01,  9.6839714164065130e-03,
    -1.0070086522629434e-01, -7.8567420131838595e-02, 4.6957584637855376e-02,
    1.1068829978797175e-01,  2.8757671678057914e-02,  -9.1016893809887853e-02,
    8.5116049235441998e-02,  -5.5555555555555532e-02, -1.0441029119843427e-01,
    1.9294241962992217e-02,  1.1111111111111110e-01,  1.9294241962992370e-02,
    -1.0441029119843427e-01, -5.5555555555555705e-02, 8.5116049235442026e-02,
    8.5116049235442096e-02,  -5.5555555555555615e-02, -1.0441029119843431e-01,
    1.9294241962992360e-02,  1.1111111111111110e-01,  1.9294241962992134e-02,
    -1.0441029119843424e-01, -5.5555555555555414e-02, 8.5116049235441985e-02,
    7.8567420131838622e-02,  -7.8567420131838608e-02, -7.8567420131838636e-02,
    7.8567420131838595e-02,  7.8567420131838636e-02,  -7.8567420131838650e-02,
    -7.8567420131838567e-02, 7.8567420131838511e-02,  7.8567420131838581e-02,
    -7.8567420131838497e-02, -7.8567420131838595e-02, 7.8567420131838483e-02,
    7.8567420131838608e-02,  -7.8567420131838483e-02, -7.8567420131838886e-02,
    7.8567420131838470e-02,  7.8567420131838622e-02,  -7.8567420131838747e-02,
    7.1420845520726597e-02,  -9.6225044864937603e-02, -3.8002238147296502e-02,
    1.0942308366802311e-01,  3.4017966642982035e-17,  -1.0942308366802313e-01,
    3.8002238147296440e-02,  9.6225044864937687e-02,  -7.1420845520726625e-02,
    -7.1420845520726736e-02, 9.6225044864937617e-02,  3.8002238147296384e-02,
    -1.0942308366802310e-01, -2.9942688208452948e-16, 1.0942308366802314e-01,
    -3.8002238147296558e-02, -9.6225044864937923e-02, 7.1420845520726278e-02,
    6.3730715150116246e-02,  -1.0732509180989647e-01, 9.6839714164064315e-03,
    1.0070086522629447e-01,  -7.8567420131838525e-02, -4.6957584637855668e-02,
    1.1068829978797171e-01,  -2.8757671678057581e-02, -9.1016893809888172e-02,
    9.1016893809887867e-02,  2.8757671678058299e-02,  -1.1068829978797175e-01,
    4.6957584637854995e-02,  7.8567420131838900e-02,  -1.0070086522629433e-01,
    -9.6839714164069744e-03, 1.0732509180989656e-01,  -6.3730715150115719e-02,
    5.5555555555555566e-02,  -1.1111111111111110e-01, 5.5555555555555483e-02,
    5.5555555555555525e-02,  -1.1111111111111110e-01, 5.5555555555555629e-02,
    5.5555555555555719e-02,  -1.1111111111111110e-01, 5.5555555555555608e-02,
    5.5555555555555740e-02,  -1.1111111111111110e-01, 5.5555555555555580e-02,
    5.5555555555555768e-02,  -1.1111111111111110e-01, 5.5555555555555219e-02,
    5.5555555555555788e-02,  -1.1111111111111110e-01, 5.5555555555555532e-02,
    4.6957584637855494e-02,  -1.0732509180989648e-01, 9.1016893809888005e-02,
    -9.6839714164066153e-03, -7.8567420131838567e-02, 1.1068829978797173e-01,
    -6.3730715150116274e-02, -2.8757671678057522e-02, 1.0070086522629441e-01,
    -1.0070086522629451e-01, 2.8757671678057730e-02,  6.3730715150115941e-02,
    -1.1068829978797172e-01, 7.8567420131838442e-02,  9.6839714164062025e-03,
    -9.1016893809887769e-02, 1.0732509180989659e-01,  -4.6957584637855293e-02,
    3.8002238147296537e-02,  -9.6225044864937617e-02, 1.0942308366802311e-01,
    -7.1420845520726500e-02, -4.7625153300174848e-17, 7.1420845520726570e-02,
    -1.0942308366802314e-01, 9.6225044864937617e-02,  -3.8002238147296218e-02,
    -3.8002238147296766e-02, 9.6225044864937520e-02,  -1.0942308366802310e-01,
    7.1420845520726278e-02,  -5.4497522255058850e-17, -7.1420845520726792e-02,
    1.0942308366802322e-01,  -9.6225044864937576e-02, 3.8002238147296120e-02,
    2.8757671678057886e-02,  -7.8567420131838692e-02, 1.0732509180989649e-01,
    -1.0732509180989641e-01, 7.8567420131838511e-02,  -2.8757671678057581e-02,
    -2.8757671678057904e-02, 7.8567420131838608e-02,  -1.0732509180989656e-01,
    1.0732509180989641e-01,  -7.8567420131838456e-02, 2.8757671678057324e-02,
    2.8757671678057966e-02,  -7.8567420131839205e-02, 1.0732509180989667e-01,
    -1.0732509180989629e-01, 7.8567420131838137e-02,  -2.8757671678057255e-02,
    1.9294241962992269e-02,  -5.5555555555555608e-02, 8.5116049235442012e-02,
    -1.0441029119843430e-01, 1.1111111111111110e-01,  -1.0441029119843420e-01,
    8.5116049235442012e-02,  -5.5555555555555247e-02, 1.9294241962992332e-02,
    1.9294241962992547e-02,  -5.5555555555555441e-02, 8.5116049235442151e-02,
    -1.0441029119843420e-01, 1.1111111111111110e-01,  -1.0441029119843435e-01,
    8.5116049235441943e-02,  -5.5555555555555837e-02, 1.9294241962992224e-02,
    9.6839714164064592e-03,  -2.8757671678057848e-02, 4.6957584637855473e-02,
    -6.3730715150116204e-02, 7.8567420131838581e-02,  -9.1016893809887950e-02,
    1.0070086522629441e-01,  -1.0732509180989645e-01, 1.1068829978797172e-01,
    -1.1068829978797173e-01, 1.0732509180989661e-01,  -1.0070086522629432e-01,
    9.1016893809888269e-02,  -7.8567420131838428e-02, 6.3730715150116024e-02,
    -4.6957584637855980e-02, 2.8757671678057626e-02,  -9.6839714164070143e-03};

// Pact storage of lower triangular transform matrix. (18 x 18 + 18) / 2 items
const double CM[171] = {
    1.0000000000000000e+00,  -5.0000000000000000e-01, 5.0000000000000000e-01,
    -2.5000000000000000e-01, -1.0000000000000000e+00, 2.5000000000000000e-01,
    2.5000000000000000e-01,  3.7500000000000000e-01,  -7.5000000000000000e-01,
    1.2500000000000000e-01,  1.8750000000000000e-01,  5.0000000000000000e-01,
    7.5000000000000000e-01,  -5.0000000000000000e-01, 6.2500000000000000e-02,
    -1.8750000000000000e-01, -3.1250000000000000e-01, 0.0000000000000000e+00,
    7.8125000000000000e-01,  -3.1250000000000000e-01, 3.1250000000000000e-02,
    -1.5625000000000000e-01, -3.7500000000000000e-01, -5.1562500000000000e-01,
    -4.3750000000000000e-01, 6.5625000000000000e-01,  -1.8750000000000000e-01,
    1.5625000000000000e-02,  1.5625000000000000e-01,  2.7343750000000000e-01,
    1.0937500000000000e-01,  -2.7343750000000000e-01, -6.5625000000000000e-01,
    4.9218750000000000e-01,  -1.0937500000000000e-01, 7.8125000000000000e-03,
    1.3671875000000000e-01,  3.1250000000000000e-01,  4.0625000000000000e-01,
    4.3750000000000000e-01,  1.0937500000000000e-01,  -6.8750000000000000e-01,
    3.4375000000000000e-01,  -6.2500000000000000e-02, 3.9062500000000000e-03,
    -1.3671875000000000e-01, -2.4609375000000000e-01, -1.4062500000000000e-01,
    9.3750000000000000e-02,  4.2187500000000000e-01,  4.2187500000000000e-01,
    -6.0937500000000000e-01, 2.2851562500000000e-01,  -3.5156250000000000e-02,
    1.9531250000000000e-03,  -1.2304687500000000e-01, -2.7343750000000000e-01,
    -3.4179687500000000e-01, -3.9062500000000000e-01, -2.7343750000000000e-01,
    1.7187500000000000e-01,  5.9082031250000000e-01,  -4.8828125000000000e-01,
    1.4648437500000000e-01,  -1.9531250000000000e-02, 9.7656250000000000e-04,
    1.2304687500000000e-01,  2.2558593750000000e-01,  1.5039062500000000e-01,
    -1.0742187500000000e-02, -2.5781250000000000e-01, -4.3505859375000000e-01,
    -1.3964843750000000e-01, 6.2841796875000000e-01,  -3.6523437500000000e-01,
    9.1308593750000000e-02,  -1.0742187500000000e-02, 4.8828125000000000e-04,
    1.1279296875000000e-01,  2.4609375000000000e-01,  2.9882812500000000e-01,
    3.4765625000000000e-01,  3.0834960937500000e-01,  6.4453125000000000e-02,
    -3.5449218750000000e-01, -3.9257812500000000e-01, 5.7861328125000000e-01,
    -2.5976562500000000e-01, 5.5664062500000000e-02,  -5.8593750000000000e-03,
    2.4414062500000000e-04,  -1.1279296875000000e-01, -2.0947265625000000e-01,
    -1.5234375000000000e-01, -3.3325195312500000e-02, 1.5551757812500000e-01,
    3.4753417968750000e-01,  3.3007812500000000e-01,  -1.2377929687500000e-01,
    -5.3955078125000000e-01, 4.8559570312500000e-01,  -1.7773437500000000e-01,
    3.3325195312500000e-02,  -3.1738281250000000e-03, 1.2207031250000000e-04,
    -1.0473632812500000e-01, -2.2558593750000000e-01, -2.6788330078125000e-01,
    -3.1274414062500000e-01, -3.0676269531250000e-01, -1.6918945312500000e-01,
    1.3629150390625000e-01,  4.1162109375000000e-01,  1.4184570312500000e-01,
    -5.8447265625000000e-01, 3.8153076171875000e-01,  -1.1791992187500000e-01,
    1.9653320312500000e-02,  -1.7089843750000000e-03, 6.1035156250000000e-05,
    1.0473632812500000e-01,  1.9638061523437500e-01,  1.5106201171875000e-01,
    5.8746337890625000e-02,  -8.9721679687500000e-02, -2.6358032226562500e-01,
    -3.4515380859375000e-01, -1.4877319335937500e-01, 3.1127929687500000e-01,
    3.6056518554687500e-01,  -5.5499267578125000e-01, 2.8518676757812500e-01,
    -7.6293945312500000e-02, 1.1444091796875000e-02,  -9.1552734375000000e-04,
    3.0517578125000000e-05,  9.8190307617187500e-02,  2.0947265625000000e-01,
    2.4438476562500000e-01,  2.8466796875000000e-01,  2.9406738281250000e-01,
    2.1533203125000000e-01,  2.6855468750000000e-03,  -2.7978515625000000e-01,
    -3.4722900390625000e-01, 1.0205078125000000e-01,  4.9633789062500000e-01,
    -4.8291015625000000e-01, 2.0495605468750000e-01,  -4.8339843750000000e-02,
    6.5917968750000000e-03,  -4.8828125000000000e-04, 1.5258789062500000e-05,
    -9.8190307617187500e-02, -1.8547058105468750e-01, -1.4837646484375000e-01,
    -7.4188232421875000e-02, 4.5654296875000000e-02,  1.9662475585937500e-01,
    3.1024169921875000e-01,  2.5628662109375000e-01,  -5.2917480468750000e-02,
    -3.8806152343750000e-01, -1.3177490234375000e-01, 5.4837036132812500e-01,
    -3.9428710937500000e-01, 1.4266967773437500e-01,  -3.0090332031250000e-02,
    3.7612915039062500e-03,  -2.5939941406250000e-04, 7.6293945312500000e-06};

// Transpose of CM
const double CMT[171] = {
    1.0000000000000000e+00,  -5.0000000000000000e-01, -2.5000000000000000e-01,
    2.5000000000000000e-01,  1.8750000000000000e-01,  -1.8750000000000000e-01,
    -1.5625000000000000e-01, 1.5625000000000000e-01,  1.3671875000000000e-01,
    -1.3671875000000000e-01, -1.2304687500000000e-01, 1.2304687500000000e-01,
    1.1279296875000000e-01,  -1.1279296875000000e-01, -1.0473632812500000e-01,
    1.0473632812500000e-01,  9.8190307617187500e-02,  -9.8190307617187500e-02,
    5.0000000000000000e-01,  -1.0000000000000000e+00, 3.7500000000000000e-01,
    5.0000000000000000e-01,  -3.1250000000000000e-01, -3.7500000000000000e-01,
    2.7343750000000000e-01,  3.1250000000000000e-01,  -2.4609375000000000e-01,
    -2.7343750000000000e-01, 2.2558593750000000e-01,  2.4609375000000000e-01,
    -2.0947265625000000e-01, -2.2558593750000000e-01, 1.9638061523437500e-01,
    2.0947265625000000e-01,  -1.8547058105468750e-01, 2.5000000000000000e-01,
    -7.5000000000000000e-01, 7.5000000000000000e-01,  0.0000000000000000e+00,
    -5.1562500000000000e-01, 1.0937500000000000e-01,  4.0625000000000000e-01,
    -1.4062500000000000e-01, -3.4179687500000000e-01, 1.5039062500000000e-01,
    2.9882812500000000e-01,  -1.5234375000000000e-01, -2.6788330078125000e-01,
    1.5106201171875000e-01,  2.4438476562500000e-01,  -1.4837646484375000e-01,
    1.2500000000000000e-01,  -5.0000000000000000e-01, 7.8125000000000000e-01,
    -4.3750000000000000e-01, -2.7343750000000000e-01, 4.3750000000000000e-01,
    9.3750000000000000e-02,  -3.9062500000000000e-01, -1.0742187500000000e-02,
    3.4765625000000000e-01,  -3.3325195312500000e-02, -3.1274414062500000e-01,
    5.8746337890625000e-02,  2.8466796875000000e-01,  -7.4188232421875000e-02,
    6.2500000000000000e-02,  -3.1250000000000000e-01, 6.5625000000000000e-01,
    -6.5625000000000000e-01, 1.0937500000000000e-01,  4.2187500000000000e-01,
    -2.7343750000000000e-01, -2.5781250000000000e-01, 3.0834960937500000e-01,
    1.5551757812500000e-01,  -3.0676269531250000e-01, -8.9721679687500000e-02,
    2.9406738281250000e-01,  4.5654296875000000e-02,  3.1250000000000000e-02,
    -1.8750000000000000e-01, 4.9218750000000000e-01,  -6.8750000000000000e-01,
    4.2187500000000000e-01,  1.7187500000000000e-01,  -4.3505859375000000e-01,
    6.4453125000000000e-02,  3.4753417968750000e-01,  -1.6918945312500000e-01,
    -2.6358032226562500e-01, 2.1533203125000000e-01,  1.9662475585937500e-01,
    1.5625000000000000e-02,  -1.0937500000000000e-01, 3.4375000000000000e-01,
    -6.0937500000000000e-01, 5.9082031250000000e-01,  -1.3964843750000000e-01,
    -3.5449218750000000e-01, 3.3007812500000000e-01,  1.3629150390625000e-01,
    -3.4515380859375000e-01, 2.6855468750000000e-03,  3.1024169921875000e-01,
    7.8125000000000000e-03,  -6.2500000000000000e-02, 2.2851562500000000e-01,
    -4.8828125000000000e-01, 6.2841796875000000e-01,  -3.9257812500000000e-01,
    -1.2377929687500000e-01, 4.1162109375000000e-01,  -1.4877319335937500e-01,
    -2.7978515625000000e-01, 2.5628662109375000e-01,  3.9062500000000000e-03,
    -3.5156250000000000e-02, 1.4648437500000000e-01,  -3.6523437500000000e-01,
    5.7861328125000000e-01,  -5.3955078125000000e-01, 1.4184570312500000e-01,
    3.1127929687500000e-01,  -3.4722900390625000e-01, -5.2917480468750000e-02,
    1.9531250000000000e-03,  -1.9531250000000000e-02, 9.1308593750000000e-02,
    -2.5976562500000000e-01, 4.8559570312500000e-01,  -5.8447265625000000e-01,
    3.6056518554687500e-01,  1.0205078125000000e-01,  -3.8806152343750000e-01,
    9.7656250000000000e-04,  -1.0742187500000000e-02, 5.5664062500000000e-02,
    -1.7773437500000000e-01, 3.8153076171875000e-01,  -5.5499267578125000e-01,
    4.9633789062500000e-01,  -1.3177490234375000e-01, 4.8828125000000000e-04,
    -5.8593750000000000e-03, 3.3325195312500000e-02,  -1.1791992187500000e-01,
    2.8518676757812500e-01,  -4.8291015625000000e-01, 5.4837036132812500e-01,
    2.4414062500000000e-04,  -3.1738281250000000e-03, 1.9653320312500000e-02,
    -7.6293945312500000e-02, 2.0495605468750000e-01,  -3.9428710937500000e-01,
    1.2207031250000000e-04,  -1.7089843750000000e-03, 1.1444091796875000e-02,
    -4.8339843750000000e-02, 1.4266967773437500e-01,  6.1035156250000000e-05,
    -9.1552734375000000e-04, 6.5917968750000000e-03,  -3.0090332031250000e-02,
    3.0517578125000000e-05,  -4.8828125000000000e-04, 3.7612915039062500e-03,
    1.5258789062500000e-05,  -2.5939941406250000e-04, 7.6293945312500000e-06};

size_t get_number_of_blocks(const size_t level) {
  return pow(2, level + 1) - 1;
}

size_t get_h(const size_t level, const size_t L) {
  return pow(2, L - level - 1);
}

size_t get_number_of_submatrices(const size_t level) {
  return 3 * get_number_of_blocks(level);
}

size_t get_total_number_of_submatrices(const size_t L) {
  return 3 * get_total_number_of_blocks(L);
}

size_t get_total_number_of_blocks(const size_t L) {
  return pow(2, L + 1) - (L + 2);
}

void get_ij(size_t *ij, const size_t level, const size_t block, const size_t s,
            const size_t L) {
  size_t h = get_h(level, L);
  ij[0] = 2 * block * s * h;
  ij[1] = ij[0] + 2 * s * h;
}

size_t direct(const double *u, double *b, direct_plan *dplan, size_t direction,
              size_t strides) {
  size_t flops = 0;
  const size_t N = dplan->N;
  for (size_t i = 0; i < N; i++)
    b[i] = 0.0;
  flops += N;

  if (direction == L2C) {
    const double *a = dplan->a;
    for (size_t n = 0; n < N; n = n + 2) {
      const double *ap = &a[n / 2];
      const double *cp = &u[n];
      const double a0 = ap[0];
      for (size_t i = 0; i < N - n; i++) {
        b[i * strides] += a0 * ap[i] * cp[i];
      }
      flops += 3 * (N - n);
    }
    b[0] /= 2;
    for (size_t i = 0; i < N; i++) {
      b[i] *= M_2_PI;
    }
    flops += N;
  } else {
    double *vn = (double *)fftw_malloc(N * sizeof(double));
    const double *an = dplan->an;
    const double *dn = dplan->dn;
    vn[0] = u[0];
    for (size_t i = 1; i < N; i++) {
      vn[i] = u[i * strides] * i;
    }
    for (size_t n = 0; n < N; n++)
      b[n * strides] = sqrt(M_PI) * an[n] * vn[n];

    for (size_t n = 2; n < N; n = n + 2) {
      const double *ap = &an[n / 2];
      const double *vp = &vn[n];
      for (size_t i = 0; i < N - n; i++)
        b[i * strides] -= dn[n / 2] * ap[i] * vp[i];
      flops += 3 * (N - n);
    }

    for (size_t i = 0; i < N; i++)
      b[i * strides] *= (i + 0.5);
    flops += N;
    fftw_free(vn);
  }
  return flops;
}

size_t directM(const double *input_array, double *output_array,
               fmm_plan *fmmplan, const size_t strides) {
  size_t s = fmmplan->s;
  size_t N = fmmplan->N;
  const double *a = fmmplan->dplan->a;
  size_t flops = 0;
  size_t h = 2 * s;
  size_t nL = N / h;

  // This is the most straightforward implementation,
  // but this is slow for large N
  /*
   for (size_t n = 0; n < s; n++) {
     const double *ap = &a[n];
     const double a0 = ap[0];
     const double *cp = &input_array[2 * n * strides];
     double *op = &output_array[0];
     if (strides == 1) {
       for (size_t i = 0; i < N - 2 * n; i++) {
         (*op++) += a0 * (*ap++) * (*cp++);
       }
     } else {
       for (size_t i = 0; i < N - 2 * n; i++) {
         (*op) += a0 * (*ap++) * (*cp);
         op += strides;
         cp += strides;
       }
     }
     flops += (N - 2 * n) * 3;
   }
   */

  // Following implementation is faster, especially so for large N.
  // First blockwise for all but last two blocks that are complicated
  // due to N and Nn

  //////////////////
  /*
  for (size_t block = 0; block < nL - 1; block++)
  {
    size_t i0 = block * h;
    double *vp = &output_array[i0];
    for (size_t i = 0; i < 2 * s; i = i + 2)
    {
      double s0 = 0.0;
      double s1 = 0.0;
      const double *ap0 = &a[0];
      const double *ap1 = &a[i0+i];
      const double *up = &input_array[i0+i];
      for (size_t j = i; j < 4 * s; j = j + 2)
      {
        s0 += (*ap0) * (*ap1++) * (*up++);
        s1 += (*ap0++) * (*ap1) * (*up++);
      }
      (*vp++) += s0;
      (*vp++) += s1;
    }
  }
  */
  /////////////////

  for (size_t block = 0; block < nL - 1; block++) {
    size_t i0 = block * h;
    for (size_t n = 0; n < 4 * s; n = n + 2) {
      const size_t n1 = n / 2;
      const double *ap = &a[i0 + n1];
      const double a0 = a[n1];
      size_t i;
      if (strides == 1) {
        double *vp = &output_array[i0];
        const double *up = &input_array[i0 + n];
        for (i = 0; i < lmin(h, 2 * h - n); i++) {
          (*vp++) += a0 * (*ap++) * (*up++);
        }
      } else {
        double *vp = &output_array[i0 * strides];
        const double *up = &input_array[(i0 + n) * strides];
        for (i = 0; i < lmin(h, 2 * h - n); i++) {
          (*vp) += a0 * (*ap++) * (*up);
          vp += strides;
          up += strides;
        }
      }
      flops += i * 3;
    }
  }

  // Last block
  size_t i0 = (nL - 1) * h;
  for (size_t n = 0; n < N - i0; n++) {
    const double *ap = &a[n + i0];
    const double a0 = a[n];
    const double *cp = &input_array[(i0 + 2 * n) * strides];
    double *op = &output_array[i0 * strides];
    if (strides == 1) {
      for (size_t i = i0; i < N - 2 * n; i++) {
        (*op++) += a0 * (*ap++) * (*cp++);
      }
    } else {
      for (size_t i = i0; i < N - 2 * n; i++) {
        (*op) += a0 * (*ap++) * (*cp);
        op += strides;
        cp += strides;
      }
    }
    flops += (N - (2 * n + i0)) * 3;
  }

  double *op = &output_array[0];
  if (strides == 1) {
    {
#pragma omp parallel for
      for (size_t i = 0; i < N; i++)
        output_array[i] *= M_2_PI;
    }
  } else {
    for (size_t i = 0; i < N; i++) {
      (*op) *= M_2_PI;
      op += strides;
    }
  }
  output_array[0] *= 0.5;
  return flops;
}

size_t directL(const double *input, double *output_array, fmm_plan *fmmplan,
               size_t strides) {
  size_t s = fmmplan->s;
  size_t N = fmmplan->N;
  const double *an = fmmplan->dplan->an;
  const double *dn = fmmplan->dplan->dn;
  size_t h = 2 * s;
  size_t nL = N / h;
  const double SPI = sqrt(M_PI);
  size_t flops = 0;
  double *op = &output_array[0];
  const double *ia = &input[0];
  const double *ap = &an[0];
  if (strides == 1) {
    for (size_t i = 0; i < N; i++)
      (*op++) += SPI * (*ap++) * (*ia++);
  } else {
    for (size_t i = 0; i < N; i++) {
      (*op) += SPI * (*ap++) * (*ia++);
      op += strides;
    }
  }
  flops += N * 3;

  // This is the most straightforward implementation
  // for (size_t n = 1; n < s; n++) {
  //  ap = &an[n];
  //  ia = &input[2 * n];
  //  op = &output_array[0];
  //  if (strides == 1) {
  //    for (size_t i = 0; i < N - 2 * n; i++) {
  //      (*op++) -= dn[n] * (*ap++) * (*ia++);
  //    }
  //  } else {
  //    for (size_t i = 0; i < N - 2 * n; i++) {
  //      (*op) -= dn[n] * (*ap++) * (*ia++);
  //      op += strides;
  //    }
  //  }
  //  flops += (N - 2 * n) * 3;
  //}

  // This implementation is faster
  // First blockwise for all but last two blocks that are complicated
  // due to N and Nn
  for (size_t block = 0; block < nL - 1; block++) {
    size_t i0 = block * h;
    for (size_t n = 2; n < 4 * s; n = n + 2) {
      const size_t n1 = n / 2;
      const double *ap = &an[i0 + n1];
      const double d0 = dn[n1];
      size_t i;
      if (strides == 1) {
        double *vp = &output_array[i0];
        const double *ia = &input[i0 + n];
        for (i = 0; i < lmin(h, 2 * h - n); i++) {
          (*vp++) -= d0 * (*ap++) * (*ia++);
        }
      } else {
        double *vp = &output_array[i0 * strides];
        const double *ia = &input[i0 + n];
        for (i = 0; i < lmin(h, 2 * h - n); i++) {
          (*vp) -= d0 * (*ap++) * (*ia++);
          vp += strides;
        }
      }
      flops += i * 3;
    }
  }

  // Last two blocks
  size_t i0 = (nL - 1) * h;
  for (size_t n = 1; n < N - i0; n++) {
    const double *ap = &an[n + i0];
    const double *ia = &input[i0 + 2 * n];
    double *op = &output_array[i0 * strides];
    if (strides == 1) {
      for (size_t i = i0; i < N - 2 * n; i++) {
        (*op++) -= dn[n] * (*ap++) * (*ia++);
      }
    } else {
      for (size_t i = i0; i < N - 2 * n; i++) {
        (*op) -= dn[n] * (*ap++) * (*ia++);
        op += strides;
      }
    }
    flops += (N - (2 * n + i0)) * 3;
  }

  // Multiply result by (x+1/2) since we have created Chebyshev
  // approximation for L(x, y)/(x+1/2)
  op = &output_array[0];
  if (strides == 1) {
    for (size_t i = 0; i < N; i++) {
      (*op++) *= (i + 0.5);
    }
  } else {
    for (size_t i = 0; i < N; i++) {
      (*op) *= (i + 0.5);
      op += strides;
    }
  }
  return flops;
}

void matvectri(const double *A, const double *x, double *b, double *w,
                const size_t m, const bool upper) {
  // compact triangular matrix A
  size_t i, j;
  if (upper == false) {
    double *zp = &w[0];
    double *zm = &w[m];
    const double *xp = &x[m];

    if (m % 2 == 0) {
      for (i = 0; i < m; i = i + 2) {
        zp[i] = x[i] + xp[i];
        zm[i] = x[i] - xp[i];
        zp[i + 1] = x[i + 1] - xp[i + 1];
        zm[i + 1] = x[i + 1] + xp[i + 1];
      }
      const double *a0 = &A[0];
      for (i = 0; i < m; i = i + 2) {
        const double *z0 = &zp[0];
        const double *z1 = &zm[0];
        const double *a1 = a0 + i + 1;
        double s0 = (*a0++) * (*z0++);
        double s1 = (*a1++) * (*z1++);
        for (size_t j = 1; j < i + 1; j++) {
          s0 += (*a0++) * (*z0++);
          s1 += (*a1++) * (*z1++);
        }
        s1 += (*a1++) * (*z1);
        a0 = a1;
        b[i] = s0;
        b[i + 1] = s1;
      }
    } else {
      for (i = 0; i < m-2; i = i + 2) {
        zp[i] = x[i] + xp[i];
        zm[i] = x[i] - xp[i];
        zp[i + 1] = x[i + 1] - xp[i + 1];
        zm[i + 1] = x[i + 1] + xp[i + 1];
      }
      zp[i] = x[i] + xp[i];
      zm[i] = x[i] - xp[i];

      const double *a0 = &A[0];
      for (i = 0; i < m-2; i = i + 2) {
        const double *a1 = a0 + i + 1;
        double s0 = (*a0++) * zp[0];
        double s1 = (*a1++) * zm[0];
        for (j = 1; j < i + 1; j++) {
          s0 += (*a0++) * zp[j];
          s1 += (*a1++) * zm[j];
        }
        s1 += (*a1++) * zm[j];
        a0 = a1;
        b[i] = s0;
        b[i + 1] = s1;
      }
      double s0 = (*a0++) * zp[0];
      for (j = 1; j < i + 1; j++) {
        s0 += (*a0++) * zp[j];
      }
      b[i] = s0;
    }
  } else {
    double *bp = &b[m];
    const double *ap = &A[0];
    for (i = 0; i < m; i++) {
      double se = 0.0;
      double so = 0.0;
      for (j = i; j < m - 1; j = j + 2) {
        se += (*ap++) * x[j];
        so += (*ap++) * x[j+1];
      }
      if ((i+m) % 2 == 1)
        se += (*ap++) * x[j];
      (*b++) += se + so;
      (*bp++) += se - so;
    }
  }
}

void vandermonde(double *T, const size_t h, const size_t N) {
  double *x = (double *)fftw_malloc(2 * h * sizeof(double));
  double *x2 = (double *)fftw_malloc(2 * h * sizeof(double));
  double *Tm = (double *)fftw_malloc(2 * h * N * sizeof(double));

  for (size_t i = 0; i < 2 * h; i++) {
    x[i] = -1 + (double)(i) / ((double)h);
    x2[i] = 2 * x[i];
    Tm[i * N] = 1;
    Tm[i * N + 1] = x[i];
  }

  for (size_t i = 0; i < 2 * h; i++) {
    for (size_t j = 2; j < N; j++) {
      Tm[i * N + j] = Tm[i * N + j - 1] * x2[i] - Tm[i * N + j - 2];
    }
  }

  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < N; j++) {
      T[i * N + j] = Tm[2 * i * N + j];             // even
      T[(h + i) * N + j] = Tm[(2 * i + 1) * N + j]; // odd
    }
  }
  fftw_free(x);
  fftw_free(x2);
  fftw_free(Tm);
}

double sum(const double *a, const size_t N) {
  double x = 0;
  for (size_t i = 0; i < N; i++)
    x += a[i];
  return x;
}

void free_direct(direct_plan *plan) {
  if (plan->a != NULL) {
    fftw_free(plan->a);
    plan->a = NULL;
  }
  if (plan->an != NULL) {
    fftw_free(plan->an);
    plan->an = NULL;
  }
  if (plan->dn != NULL) {
    fftw_free(plan->dn);
    plan->dn = NULL;
  }
  fftw_free(plan);
  plan = NULL;
}

void free_fmm_2d(fmm_plan_2d *plan) {
  if (plan->fmmplan0 == plan->fmmplan1) {
    free_fmm(plan->fmmplan0);
    plan->fmmplan1 = NULL;
  } else if (plan->fmmplan0 != NULL) {
    free_fmm(plan->fmmplan0);
  } else if (plan->fmmplan1 != NULL) {
    free_fmm(plan->fmmplan1);
  }
  fftw_free(plan);
  plan = NULL;
}

void free_fmm(fmm_plan *plan) {
  if (plan->A[0] != NULL) {
    fftw_free(plan->A[0]);
    plan->A[0] = NULL;
  }
  if (plan->A[1] != NULL) {
    fftw_free(plan->A[1]);
    plan->A[1] = NULL;
  }
  if (plan->A != NULL) {
    fftw_free(plan->A);
    plan->A = NULL;
  }
  if (plan->T != NULL) {
    fftw_free(plan->T);
    plan->T = NULL;
  }
  if (plan->M != 18) {
    if (plan->Th != NULL) {
      fftw_free(plan->Th);
      plan->Th = NULL;
    }
    if (plan->ThT != NULL) {
      fftw_free(plan->ThT);
      plan->ThT = NULL;
    }
  }
  if (plan->ia != NULL) {
    fftw_free(plan->ia);
    plan->ia = NULL;
  }
  if (plan->oa != NULL) {
    fftw_free(plan->ia);
    plan->ia = NULL;
  }
  if (plan->work != NULL) {
    fftw_free(plan->work);
    plan->work = NULL;
  }
  if (plan->wk != NULL) {
    fftw_free(plan->wk[0]);
    for (size_t level = 0; level < plan->L; level++)
      plan->wk[level] = NULL;
    fftw_free(plan->wk);
  }
  if (plan->ck != NULL) {
    fftw_free(plan->ck[0]);
    for (size_t level = 0; level < plan->L; level++)
      plan->ck[level] = NULL;
    fftw_free(plan->ck);
  }
  if (plan->dplan != NULL) {
    free_direct(plan->dplan);
  }
  fftw_free(plan);
  plan = NULL;
}

direct_plan *create_direct(size_t N, size_t direction) {
  direct_plan *dplan = (direct_plan *)fftw_malloc(sizeof(direct_plan));
  dplan->a = NULL;
  dplan->an = NULL;
  dplan->dn = NULL;
  if (direction == L2C | direction == BOTH) {
    double *a = (double *)fftw_malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++)
      a[i] = _LambdaI(i);
    dplan->a = a;
  }
  if (direction == C2L | direction == BOTH) {
    double *dn = (double *)fftw_malloc((N + 1) / 2 * sizeof(double));
    double *an = (double *)fftw_malloc(N * sizeof(double));
    dn[0] = 0;
    an[0] = M_2_SQRTPI;
    for (size_t i = 1; i < N; i++)
      an[i] = 1 / (2 * _LambdaI(i) * i * (i + 0.5));
    for (size_t i = 1; i < (N + 1) / 2; i++)
      dn[i] = _LambdaI(i - 1) / (2 * i);

    dplan->an = an;
    dplan->dn = dn;
  }
  dplan->direction = direction;
  dplan->N = N;
  return dplan;
}

// Direct matrix-vector DCT of fixed size N=18 and precomputed matrix is
// faster than FFTW
void dct(double *input, double *output) {
  const size_t N = 18;
  cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1, &DCTII[0], N, &input[0], 1,
              0, &output[0], 1);
}

void dct2(double *input, double *output) {
  const size_t N = 18;
  double *tmp = (double *)fftw_malloc(N * N * sizeof(double));

  // DCT along rows
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, &input[0],
              N, &DCTII[0], N, 0, &tmp[0], N);

  // DCT along columns
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0,
              &DCTII[0], N, &tmp[0], N, 0, &output[0], N);

  fftw_free(tmp);
}

fmm_plan *create_fmm(const size_t N, const size_t maxs, const size_t M,
                     const size_t direction, size_t lagrange, const size_t v) {
  fmm_plan *fmmplan = (fmm_plan *)fftw_malloc(sizeof(fmm_plan));
  fftw_plan plan1d, plan;
  size_t Nn;
  size_t s;
  size_t ij[2];
  size_t directions[2];
  size_t num_directions = 2;
  assert(maxs >= 32); // Because of implementation for _Lambda/_LambdaE. Easily adjusted for less, but there should be no need.
  switch (direction) {
  case L2C:
    directions[0] = 0;
    num_directions = 1;
    break;
  case C2L:
    directions[0] = 1;
    num_directions = 1;
    break;
  default:
    directions[0] = 0;
    directions[1] = 1;
    break;
  }
  double **A = (double **)calloc(2, sizeof(double *));
  fmmplan->A = A;
  fmmplan->T = NULL;
  fmmplan->Th = NULL;
  fmmplan->ThT = NULL;
  fmmplan->ia = NULL;
  fmmplan->oa = NULL;
  fmmplan->work = NULL;
  fmmplan->wk = NULL;
  fmmplan->ck = NULL;
  fmmplan->dplan = NULL;
  fmmplan->lagrange = lagrange;

  int L = ceil(log2((double)N / (double)maxs)) - 2;
  if (L < 1) {
    if (v > 1)
      printf("Levels < 1. Using only direct method\n");
    fmmplan->dplan = create_direct(N, direction);
    fmmplan->Nn = N;
    fmmplan->L = 0;
    return fmmplan;
  }

  s = ceil((double)N / (double)pow(2, L + 2));
  Nn = s * pow(2, L + 2);
  fmmplan->dplan = create_direct(N, direction);
  fmmplan->M = M;
  fmmplan->L = L;
  fmmplan->N = N;
  fmmplan->Nn = Nn;
  fmmplan->s = s;
  if (v > 1) {
    printf("N %lu\n", N);
    printf("Num levels %d\n", L);
    printf("Num submatrices %lu\n", get_total_number_of_submatrices(L));
    printf("Num blocks %lu\n", get_total_number_of_blocks(L));
    printf("Given max s %lu \n", maxs);
    printf("Computed s %lu \n", s);
    printf("Computed N %lu\n", Nn);
  }

  const size_t MM = M * M;
  uint64_t t1 = tic;
  if (direction == BOTH) {
    A[0] = (double *)fftw_malloc(get_total_number_of_submatrices(L) * MM *
                                 sizeof(double));
    A[1] = (double *)fftw_malloc(get_total_number_of_submatrices(L) * MM *
                                 sizeof(double));
  } else {
    A[direction] = (double *)fftw_malloc(get_total_number_of_submatrices(L) *
                                         MM * sizeof(double));
  }

  double *xj = (double *)fftw_malloc(M * sizeof(double));
  for (size_t i = 0; i < M; i++)
    xj[i] = cos((i + 0.5) * M_PI / M);

  double *fun = (double *)fftw_malloc(MM * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(MM * sizeof(double));
  double *fx0 = (double *)fftw_malloc(2 * MM * sizeof(double));
  double *fx1 = (double *)fftw_malloc(2 * MM * sizeof(double));
  double *lx1 = (double *)fftw_malloc(2 * MM * sizeof(double));
  double *ap0 = A[0];
  double *ap1 = A[1];
  bool use_FFTW = (M != 18 && lagrange == 0);

  if (use_FFTW) {
    if (v > 1)
      printf("using FFTW for planning\n");
    plan1d = fftw_plan_r2r_1d(M, fun, fun_hat, FFTW_REDFT10, FFTW_MEASURE);
    plan = fftw_plan_r2r_2d(M, M, fun, fun_hat, FFTW_REDFT10, FFTW_REDFT10,
                            FFTW_MEASURE);
  }

  for (size_t level = 0; level < L; level++) {
    size_t h = s * get_h(level, L);
    for (size_t block = 0; block < get_number_of_blocks(level); block++) {
      get_ij(ij, level, block, s, L);
      for (size_t q = 0; q < 2; q++) {
        double y0 = 2 * (ij[1] + q * h) + h;
        for (size_t p = 0; p < q + 1; p++) {
          double x0 = 2 * (ij[0] + p * h) + h;
          for (size_t di = 0; di < num_directions; di++) {
            size_t dir = directions[di];
            for (size_t i = 0; i < M; i++) {
              double x = x0 + xj[i] * h;
              for (size_t j = 0; j < M; j++) {
                double y = y0 + xj[j] * h;
                size_t ix = q * MM + i * M + j;
                size_t xi = q * MM + j * M + i;
                size_t qp = (q - p) * MM + i * M + j;
                if (di == 0) {
                  if (block == 0 && p == 0) {
                    if (j < M - i) { // Persymmetric Lambda((y-x)/2)
                      //double m0 = _Lambda((y - x) / 2);
                      double m0 = _LambdaE((y - x) / 2);
                      fx0[ix] = m0;
                      fx0[q * MM + (M - j) * M - i - 1] = m0;
                    }
                  }
                }
                if (dir == L2C) {
                  if (j >= i) { // Symmetric Lambda((x+y)/2)
                    //double m1 = _Lambda((x + y) / 2);
                    double m1 = _LambdaE((x + y) / 2);
                    fx1[ix] = m1;
                    fx1[xi] = m1;
                  }
                  fun[i * M + j] = fx0[qp] * fx1[ix];
                } else {
                  if (j >= i) {
                    //double l1 = -_Lambda((x + y - 1) / 2) / (y + x + 1);
                    double l1 = -_LambdaE((x + y - 1) / 2) / (y + x + 1);
                    lx1[ix] = l1;
                    lx1[xi] = l1;
                  }
                  fun[i * M + j] = fx0[qp] * lx1[ix] / (y - x - 1);
                }
              }
            }
            double *fh = NULL;
            if (lagrange == 0) {
              if (use_FFTW) {
                fftw_execute(plan);
                for (size_t i = 0; i < M; i++) {
                  for (size_t j = 0; j < M; j++) {
                    fun_hat[i * M + j] /= MM;
                  }
                }
                for (size_t j = 0; j < M; j++)
                  fun_hat[j] /= 2;
                for (size_t i = 0; i < M; i++)
                  fun_hat[i * M] /= 2;
              } else {
                dct2(&fun[0], &fun_hat[0]);
              }
              fh = &fun_hat[0];
            } else {
              fh = &fun[0];
            }

            if (dir == L2C) {
              for (size_t i = 0; i < MM; i++)
                *ap0++ = *fh++;
            } else {
              for (size_t i = 0; i < MM; i++)
                *ap1++ = *fh++;
            }
          }
        }
      }
    }
  }
  double *wj = NULL;
  if (lagrange == 1) {
    wj = (double *)fftw_malloc(M * sizeof(double));
    for (size_t j = 0; j < M; j++) {
      int sign = (j % 2 == 0) ? 1 : -1;
      wj[j] = sign * sin((j + 0.5) * M_PI / M);
    }
  }

  double *T = (double *)fftw_malloc(2 * s * M * sizeof(double));
  if (lagrange == 0)
    vandermonde(T, s, M);
  else {
    double *xh = (double *)fftw_malloc(2 * s * sizeof(double));
    for (size_t i = 0; i < 2 * s; i++) {
      xh[i] = -1 + (double)(i) / ((double)s);
    }
    for (size_t i = 0; i < s; i++) {
      double sume = 0.0;
      double sumo = 0.0;
      for (size_t j = 0; j < M; j++) {
        double se = wj[j] / (xh[2 * i] - xj[j]);
        double so = wj[j] / (xh[2 * i + 1] - xj[j]);
        T[i * M + j] = se;       // even
        T[(s + i) * M + j] = so; // odd
        sume += se;
        sumo += so;
      }
      for (size_t j = 0; j < M; j++) {
        T[i * M + j] /= sume;
        T[(s + i) * M + j] /= sumo;
      }
    }
    fftw_free(xh);
  }
  fmmplan->T = T;

  double *Th = NULL;
  double *ThT = NULL;
  if (!use_FFTW && lagrange == 0) {
    if (v > 1)
      printf("Using exact binomial matrix\n");
    Th = (double *)&CM[0];
    ThT = (double *)&CMT[0];
  } else {
    if (lagrange == 0) {
      Th = (double *)fftw_malloc((MM + M) / 2 * sizeof(double));
      ThT = (double *)fftw_malloc((MM + M) / 2 * sizeof(double));
      double *Tha = (double *)fftw_malloc(MM * sizeof(double));
      double *ThTa = (double *)fftw_malloc(MM * sizeof(double));
      double *th = &Tha[0];
      for (size_t k = 0; k < M; k++) {
        for (size_t j = 0; j < M; j++) {
          fun[j] = cos(k * acos((xj[j] - 1) / 2));
          fun_hat[j] = 0.0;
        }
        fftw_execute(plan1d);
        *th++ = fun_hat[0] / M / 2;
        for (size_t j = 1; j < M; j++)
          *th++ = fun_hat[j] / M;
      }
      // Make transpose
      th = &Tha[0];
      double *tht = &ThTa[0];
      for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
          tht[i * M + j] = th[j * M + i];
        }
      }
      /// Move to compact storage
      th = &Th[0];
      for (size_t k = 0; k < M; k++) {
        for (size_t j = 0; j < k + 1; j++) {
          (*th++) = Tha[k * M + j];
        }
      }
      tht = &ThT[0];
      for (size_t k = 0; k < M; k++) {
        for (size_t j = k; j < M; j++) {
          (*tht++) = ThTa[k * M + j];
        }
      }
      fftw_free(Tha);
      fftw_free(ThTa);

    } else {
      Th = (double *)fftw_malloc(2 * MM * sizeof(double));
      double *xh = (double *)fftw_malloc(2 * M * sizeof(double));
      for (size_t i = 0; i < M; i++) {
        xh[i] = (xj[i] - 1) / 2;
        xh[M + i] = (xj[i] + 1) / 2;
      }
      for (size_t i = 0; i < M; i++) {
        double sum0 = 0.0;
        double sum1 = 0.0;
        for (size_t j = 0; j < M; j++) {
          double s0 = wj[j] / (xh[i] - xj[j]);
          double s1 = wj[j] / (xh[M + i] - xj[j]);
          Th[i * M + j] = s0;
          Th[MM + i * M + j] = s1;
          sum0 += s0;
          sum1 += s1;
        }
        for (size_t j = 0; j < M; j++) {
          Th[i * M + j] /= sum0;
          Th[MM + i * M + j] /= sum1;
        }
      }
      fftw_free(xh);
    }
  }

  double *ia = (double *)fftw_malloc(Nn / 2 * sizeof(double));
  double *oa = (double *)fftw_malloc(Nn / 2 * sizeof(double));
  double *work = (double *)fftw_malloc(2 * M * sizeof(double));
  double **wk = (double **)fftw_malloc(L * sizeof(double *));
  double **ck = (double **)fftw_malloc(L * sizeof(double *));
  size_t Nb = get_total_number_of_blocks(L);
  wk[0] = (double *)fftw_malloc(Nb * 2 * M * sizeof(double));
  ck[0] = (double *)fftw_malloc(Nb * 2 * M * sizeof(double));
  for (size_t level = 1; level < L; level++) {
    size_t b = get_number_of_blocks(level - 1);
    wk[level] = wk[level - 1] + b * 2 * M;
    ck[level] = ck[level - 1] + b * 2 * M;
  }
  fmmplan->ia = ia;
  fmmplan->oa = oa;
  fmmplan->wk = wk;
  fmmplan->ck = ck;
  fmmplan->work = work;
  if (use_FFTW) {
    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan1d);
  }

  fftw_free(fun);
  fftw_free(fun_hat);
  fftw_free(fx0);
  fftw_free(fx1);
  fftw_free(lx1);
  fftw_free(xj);
  if (lagrange == 1) {
    fftw_free(wj);
  }
  fmmplan->Th = Th;
  fmmplan->ThT = ThT;
  uint64_t t2 = tic;
  if (v > 1)
    printf("Initialization %2.4e s\n", dtics(t1, t2));
  return fmmplan;
}

fmm_plan_2d *create_fmm_2d(size_t N0, size_t N1, int axis, size_t maxs,
                           size_t M, size_t direction, size_t lagrange,
                           size_t v) {
  fmm_plan_2d *fmmplan2d = (fmm_plan_2d *)fftw_malloc(sizeof(fmm_plan_2d));
  fmmplan2d->fmmplan0 = NULL;
  fmmplan2d->fmmplan1 = NULL;
  if (v > 1) {
    printf("crate_fmm_2d\n");
  }
  if (axis == 0) {
    fmmplan2d->fmmplan0 = create_fmm(N0, maxs, M, direction, lagrange, v);
  } else if (axis == 1) {
    fmmplan2d->fmmplan1 = create_fmm(N1, maxs, M, direction, lagrange, v);
  } else if (axis == -1) {
    fmmplan2d->fmmplan0 = create_fmm(N0, maxs, M, direction, lagrange, v);
    if (N0 == N1) {
      fmmplan2d->fmmplan1 = fmmplan2d->fmmplan0;
    } else {
      fmmplan2d->fmmplan1 = create_fmm(N1, maxs, M, direction, lagrange, v);
    }
  }
  fmmplan2d->N0 = N0;
  fmmplan2d->N1 = N1;
  fmmplan2d->axis = axis;
  return fmmplan2d;
}

size_t execute2D(const double *input_array, double *output_array,
                 fmm_plan_2d *fmmplan2d, size_t direction) {
  size_t flops = 0;
  if (fmmplan2d->axis == 0) {
    for (size_t i = 0; i < fmmplan2d->N1; i++) {
      flops += execute(&input_array[i], &output_array[i], fmmplan2d->fmmplan0,
                       direction, fmmplan2d->N1);
    }
  } else if (fmmplan2d->axis == 1) {
    for (size_t i = 0; i < fmmplan2d->N0; i++) {
      size_t N1 = fmmplan2d->N1;
      flops += execute(&input_array[i * N1], &output_array[i * N1],
                       fmmplan2d->fmmplan1, direction, 1);
    }
  } else if (fmmplan2d->axis == -1) {
    double *out =
        (double *)calloc(fmmplan2d->N0 * fmmplan2d->N1, sizeof(double));
    for (size_t i = 0; i < fmmplan2d->N1; i++) {
      flops += execute(&input_array[i], &out[i], fmmplan2d->fmmplan0, direction,
                       fmmplan2d->N1);
    }

    for (size_t i = 0; i < fmmplan2d->N0; i++) {
      size_t N1 = fmmplan2d->N1;
      flops += execute(&out[i * N1], &output_array[i * N1], fmmplan2d->fmmplan1,
                       direction, 1);
    }
    fftw_free(out);
  }
  return flops;
}

size_t execute(const double *input_array, double *output_array,
               fmm_plan *fmmplan, size_t direction, const size_t stride) {
  size_t Nn = fmmplan->Nn;
  size_t N = fmmplan->N;
  size_t L = fmmplan->L;
  size_t s = fmmplan->s;
  size_t M = fmmplan->M;
  double *T = fmmplan->T;
  double *Th = fmmplan->Th;
  double *ThT = fmmplan->ThT;
  double *A = fmmplan->A[direction];
  size_t flops = 0;
  size_t lagrange = fmmplan->lagrange;
  assert(direction == C2L | direction == L2C);

  if (T == NULL) {
    flops =
        direct(input_array, output_array, fmmplan->dplan, direction, stride);
    return flops;
  }

  double *ia = fmmplan->ia;
  double *oa = fmmplan->oa;
  double **wk = fmmplan->wk;
  double **ck = fmmplan->ck;
  double *input = NULL;
  if (direction == C2L) { // Need to modify input array, so make copy
    input = (double *)fftw_malloc(N * sizeof(double));
    input[0] = input_array[0];
    input[1] = input_array[stride];
    double *w0 = &input[2];
    size_t ii = 2 * stride;
    for (size_t i = 2; i < N; i++) {
      (*w0++) = input_array[ii] * i;
      ii += stride;
    }
  }

  for (size_t odd = 0; odd < 2; odd++) {
    for (size_t i = 0; i < Nn / 2; i++) {
      oa[i] = 0.0;
    }

    double *w0 = &wk[0][0];
    double *c0 = &ck[0][0];
    for (size_t i = 0; i < 2 * M * get_total_number_of_blocks(L); i++) {
      (*w0++) = 0.0;
      (*c0++) = 0.0;
    }

    const double *ap;
    switch (direction) {
    case L2C:
      ap = &input_array[odd * stride];
      break;

    case C2L:
      ap = &input[odd];
      break;
    }

    size_t rest = N % 2;
    double *iap = &ia[0];
    if (stride == 1 || direction == C2L) {
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *iap++ = *ap;
        ap += 2;
      }
    } else {
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *iap++ = *ap;
        ap += 2 * stride;
      }
    }

    for (size_t i = N / 2 + rest * (1 - odd); i < Nn / 2; i++) {
      *iap++ = 0;
    }

    const size_t MM = M * M;
    const size_t K = get_number_of_blocks(L - 1) * 2;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, M, s, 1.0,
                &ia[2 * s], s, &T[odd * s * M], M, 0, wk[L - 1], M);
    flops += 2 * K * s * M;
    for (size_t level = L; level-- > 1;) {
      double *w1 = wk[level - 1];
      for (size_t block = 1; block < get_number_of_blocks(level); block++) {
        size_t Nd = block * 2 * M;
        double *wq = &wk[level][Nd];
        int b0 = (block - 1) / 2;
        int q0 = (block - 1) % 2;
        if (lagrange == 0) {
          matvectri(&Th[0], wq, &w1[(b0 * 2 + q0) * M], fmmplan->work, M,
                     false);
          flops += MM; //+2*M;
        } else {
          cblas_dgemv(CblasRowMajor, CblasTrans, M, M, 1, &Th[0], M, &wq[0], 1,
                      0, &w1[(b0 * 2 + q0) * M], 1);
          cblas_dgemv(CblasRowMajor, CblasTrans, M, M, 1, &Th[MM], M, &wq[M], 1,
                      0, &w1[(b0 * 2 + q0) * M], 1);

          flops += 4 * MM;
        }
      }
    }

    size_t ik = 0;
    for (size_t level = 0; level < L; level++) {
      for (size_t block = 0; block < get_number_of_blocks(level); block++) {
        size_t Nd = block * 2 * M;
        double *cp = &ck[level][Nd];
        double *wq = &wk[level][Nd];
        /*for (size_t q = 0; q < 2; q++) {
          for (size_t p = 0; p < q + 1; p++) {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &A[ik * MM],
                        M, &wq[q * M], 1, 1, &cp[p * M], 1);
            flops += 2 * MM;
            ik++;
          }
        }*/
        cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &A[ik * MM], M, wq, 1,
                    0, cp, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &A[(ik + 1) * MM], M,
                    &wq[M], 1, 1, cp, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &A[(ik + 2) * MM], M,
                    &wq[M], 1, 0, &cp[M], 1);
        flops += 6 * MM;
        ik += 3;
      }
    }

    for (size_t level = 0; level < L - 1; level++) {
      double *c0 = ck[level];
      double *c1 = ck[level + 1];
      for (size_t block = 0; block < get_number_of_blocks(level + 1) - 1;
           block++) {
        if (lagrange == 0) {
          matvectri(&ThT[0], &c0[block * M], &c1[block * 2 * M], NULL, M,
                     true);
          flops += MM;
        } else {
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &Th[0], M,
                      &c0[block * M], 1, 1, &c1[block * 2 * M], 1);
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &Th[MM], M,
                      &c0[block * M], 1, 1, &c1[block * 2 * M + M], 1);
          flops += 4 * MM;
        }
      }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, K, s, M, 1.0,
                ck[L - 1], M, &T[odd * s * M], M, 0, &oa[0], s);

    flops += 2 * K * s * M;
    double *oaa = &output_array[odd * stride];
    double *oap = &oa[0];
    if (stride == 1) {
      for (size_t i = 0; i < N; i = i + 2) {
        *oaa += (*oap++);
        oaa += 2;
      }
    } else {
      const size_t s2 = 2 * stride;
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *oaa += (*oap++);
        oaa += s2;
      }
    }
    // flops += 3*N/2;
  }

  switch (direction) {
  case L2C:
    flops += directM(input_array, output_array, fmmplan, stride);
    break;

  case C2L:
    flops += directL(input, output_array, fmmplan, stride);
    break;
  }

  if (input != NULL)
    fftw_free(input);
  return flops;
}
