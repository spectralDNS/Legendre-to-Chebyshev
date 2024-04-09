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

inline static double _LambdaE16(const double z) {
  const static double a0[5] = {9.9996949706645188e-01, -3.0498055973040047e-05,
                               4.8752082058744936e-09, -2.3645121623003757e-12,
                               2.4013727207097398e-15};
  const double z0 = z + 0.25;
  return chebval(-1 + 512 / pow(z0, 2), a0, 5) / sqrt(z0);
}

inline static double _LambdaE32(const double z) {
  const static double a0[4] = {9.9999237152186726e-01, -7.6281727285568499e-06,
                               3.0536699868997274e-10, -3.7171802006655959e-14};
  const double z0 = z + 0.25;
  return chebval(-1 + 2048 / pow(z0, 2), a0, 4) / sqrt(z0);
}

inline static double _LambdaE128(const double z) {
  const static double a0[3] = {9.9999952316642282e-01, -4.7683238349233873e-07,
                               1.1936572376807404e-12};
  const double z0 = z + 0.25;
  return chebval(-1 + 32768 / pow(z0, 2), a0, 3) / sqrt(z0);
}

inline static double _LambdaE1400(const double z) {
  const static double a0[2] = {9.9999999601403089e-01, -3.9859690541081902e-09};
  const double z0 = z + 0.25;
  // return chebval(-1 + 3920000 / pow(z0, 2), a0, 2) / sqrt(z0);
  return (a0[0] + a0[1] * (-1 + 3920000 / pow(z0, 2))) / sqrt(z0);
}

double _LambdaE(const double z) {
  if (z > 1400) {
    return _LambdaE1400(z);
  } else if (z > 128) {
    return _LambdaE128(z);
  } else if (z > 32) {
    return _LambdaE32(z);
  } else if (z > 16) {
    return _LambdaE16(z);
  }
  return _Lambda0(z);
}

double _Lambda(const double z) {
  double zz, y;
  const double zp = z + 0.25;
  const double s = 1 / sqrt(zp);
  const double z2 = 1 / (zp * zp);

  y = 1 - 1.5625e-02 * z2;
  if (z > 2014) // 2000
    goto scalereturn;
  zz = z2 * z2;
  y += 2.5634765625e-03 * zz;
  if (z > 141) // 92
    goto scalereturn;
  zz *= z2;
  y -= 1.2798309326171875e-03 * zz;
  if (z > 40) // 31
    goto scalereturn;
  zz *= z2;
  y += 1.3435110449790955e-03 * zz;
  if (z > 20)
    goto scalereturn;
  zz *= z2;
  y -= 2.4328966392204165e-03 * zz;
  if (z > 14)
    goto scalereturn;
  zz *= z2;
  y += 6.754237533641572e-03 * zz;
  if (z > 10)
    goto scalereturn;
  return _Lambda0(z);
scalereturn:
  return y * s;
}

// Lambda with integer argument. For initializing direct solver
double _LambdaI(const size_t z) {
  const static double LamInt[256] = {
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
      1.2573843161391765e-01, 1.2475610011693392e-01, 1.2379643780834211e-01,
      1.2285858600676376e-01, 1.2194173088731031e-01, 1.2104510051313890e-01,
      1.2016796210362340e-01, 1.1930961951716895e-01, 1.1846941092901987e-01,
      1.1764670668645723e-01, 1.1684090732559109e-01, 1.1605144173555330e-01,
      1.1527776545731629e-01, 1.1451935910562341e-01, 1.1377572690363885e-01,
      1.1304639532092321e-01, 1.1233091180623384e-01, 1.1162884360744486e-01,
      1.1093977667159645e-01, 1.1026331461872085e-01, 1.0959907778366831e-01,
      1.0894670232067029e-01, 1.0830583936584282e-01, 1.0767615425325071e-01,
      1.0705732578053088e-01, 1.0644904552041423e-01, 1.0585101717479392e-01,
      1.0526295596826729e-01, 1.0468458807833175e-01, 1.0411565009964517e-01,
      1.0355588853996966e-01, 1.0300505934560812e-01, 1.0246292745431544e-01,
      1.0192926637382421e-01, 1.0140385778426843e-01, 1.0088649116292012e-01,
      1.0037696342977405e-01, 9.9875078612625179e-02, 9.9380647530384461e-02,
      9.8893487493470808e-02, 9.8413422020201535e-02, 9.7940280568181340e-02,
      9.7473898279761426e-02, 9.7014115740705953e-02, 9.6560778751263399e-02,
      9.6113738108896438e-02, 9.5672849401974888e-02, 9.5237972813784100e-02,
      9.4808972936244532e-02, 9.4385718592779153e-02, 9.3968082669802250e-02,
      9.3555941956338207e-02, 9.3149176991310659e-02, 9.2747671918072247e-02,
      9.2351314345772789e-02, 9.1959995217189006e-02, 9.1573608682663010e-02,
      9.1192051979818570e-02, 9.0815225318744947e-02, 9.0443031772356644e-02,
      9.0075377171656007e-02, 8.9712170005641273e-02, 8.9353321325618698e-02,
      8.8998744653691647e-02, 8.8648355895212541e-02, 8.8302073254996866e-02,
      8.7959817157109280e-02, 8.7621510168043482e-02, 8.7287076923127288e-02,
      8.6956444055994217e-02, 8.6629540130971683e-02, 8.6306295578244180e-02,
      8.5986642631658089e-02, 8.5670515269041708e-02, 8.5357849154921117e-02,
      8.5048581585519228e-02, 8.4742651435931030e-02, 8.4439999109374136e-02,
      8.4140566488418903e-02, 8.3844296888107572e-02, 8.3551135010876423e-02,
      8.3261026903199767e-02, 8.2973919913878397e-02, 8.2689762653899351e-02,
      8.2408504957797654e-02, 8.2130097846453740e-02, 8.1854493491264307e-02,
      8.1581645179626752e-02, 8.1311507281680975e-02, 8.1044035218254387e-02,
      8.0779185429959446e-02, 8.0516915347394635e-02, 8.0257183362403048e-02,
      7.9999948800344056e-02, 7.9745171893336589e-02, 7.9492813754433622e-02,
      7.9242836352690124e-02, 7.8995202489087965e-02, 7.8749875773283351e-02,
      7.8506820601143584e-02, 7.8266002133041912e-02, 7.8027386272880209e-02,
      7.7790939647810864e-02, 7.7556629588630716e-02, 7.7324424110820439e-02,
      7.7094291896204911e-02, 7.6866202275210224e-02, 7.6640125209694890e-02,
      7.6416031276333216e-02, 7.6193891650529921e-02, 7.5973678090846306e-02,
      7.5755362923918587e-02, 7.5538919029850243e-02, 7.5324319828060898e-02,
      7.5111539263574847e-02, 7.4900551793733353e-02, 7.4691332375315098e-02,
      7.4483856452050343e-02, 7.4278099942514289e-02, 7.4074039228386498e-02,
      7.3871651143063044e-02, 7.3670912960609070e-02, 7.3471802385039850e-02,
      7.3274297539918778e-02, 7.3078376958261235e-02, 7.2884019572733952e-02,
      7.2691204706139420e-02, 7.2499912062175889e-02, 7.2310121716463394e-02,
      7.2121814107826768e-02, 7.1934970029827211e-02, 7.1749570622533843e-02,
      7.1565597364527347e-02, 7.1383032065128041e-02, 7.1201856856840912e-02,
      7.1022054188010511e-02, 7.0843606815678820e-02, 7.0666497798639621e-02,
      7.0490710490682812e-02, 7.0316228534022709e-02, 7.0143035852904420e-02,
      6.9971116647382606e-02, 6.9800455387267035e-02, 6.9631036806229979e-02,
      6.9462845896070005e-02, 6.9295867901127531e-02, 6.9130088312847310e-02,
      6.8965492864483391e-02, 6.8802067525941965e-02, 6.8639798498758134e-02,
      6.8478672211202365e-02, 6.8318675313512642e-02, 6.8159794673248661e-02,
      6.8002017370764292e-02, 6.7845330694794787e-02, 6.7689722138155342e-02,
      6.7535179393547681e-02, 6.7381690349471446e-02, 6.7229243086237345e-02,
      6.7077825872079153e-02, 6.6927427159361480e-02, 6.6778035580880774e-02,
      6.6629639946256591e-02, 6.6482229238410892e-02, 6.6335792610132449e-02,
      6.6190319380724269e-02, 6.6045799032731417e-02, 6.5902221208747211e-02,
      6.5759575708295381e-02, 6.5617852484786132e-02, 6.5477041642544101e-02,
      6.5337133433906180e-02, 6.5198118256387230e-02, 6.5059986649911833e-02,
      6.4922729294110332e-02, 6.4786337005677333e-02, 6.4650800735790978e-02,
      6.4516111567591405e-02, 6.4382260713716735e-02, 6.4249239513895010e-02,
      6.4117039432590700e-02, 6.3985652056704242e-02, 6.3855069093323211e-02,
      6.3725282367523783e-02, 6.3596283820221103e-02, 6.3468065506067428e-02,
      6.3340619591396613e-02, 6.3213938352213811e-02, 6.3088014172229326e-02,
      6.2962839540935220e-02, 6.2838407051723888e-02, 6.2714709400047267e-02,
      6.2591739381615802e-02};
  if (z > 255)
    return _Lambda((const double)z);
  return LamInt[z];
}

// Pact storage of lower triangular transform matrix. (18 x 18 + 18) / 2 items
const double BMe[171] __attribute__((aligned)) = {
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

// Transpose of BMe
const double BMeT[171] __attribute__((aligned)) = {
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

void dct2(double *input, double *output) {
  double tmp[324];
  for (size_t i = 0; i < 18; i++) {
    dct(input + i * 18, tmp + i * 18, 1);
  }
  for (size_t i = 0; i < 18; i++) {
    dct(tmp + i, output + i, 18);
  }
}

// Hard-coded DCT II for N=18
void dct(double *input, double *output, size_t st) {
  double out[18], z[9], zpm[9];
  static const double DCTIIHH0[25] __attribute__((aligned)) = {
      1.1111111111111110e-01,  1.1111111111111110e-01,  1.1111111111111110e-01,
      1.1111111111111110e-01,  1.1111111111111110e-01,  1.0441029119843427e-01,
      5.5555555555555552e-02,  -1.9294241962992262e-02, -8.5116049235441998e-02,
      -1.1111111111111110e-01, 8.5116049235441998e-02,  -5.5555555555555552e-02,
      -1.0441029119843427e-01, 1.9294241962992262e-02,  1.1111111111111110e-01,
      5.5555555555555552e-02,  -1.1111111111111110e-01, 5.5555555555555552e-02,
      5.5555555555555552e-02,  -1.1111111111111110e-01, 1.9294241962992262e-02,
      -5.5555555555555552e-02, 8.5116049235441998e-02,  -1.0441029119843427e-01,
      1.1111111111111110e-01};
  static const double DCTIIHH1[16] __attribute__((aligned)) = {
      1.0942308366802311e-01,  9.6225044864937631e-02,  7.1420845520726597e-02,
      3.8002238147296523e-02,  9.6225044864937631e-02,  0.0000000000000000e+00,
      -9.6225044864937631e-02, -9.6225044864937631e-02, 7.1420845520726597e-02,
      -9.6225044864937631e-02, -3.8002238147296523e-02, 1.0942308366802311e-01,
      3.8002238147296523e-02,  -9.6225044864937631e-02, 1.0942308366802311e-01,
      -7.1420845520726597e-02};
  static const double DCTIIa[16] __attribute__((aligned)) = {
      1.1068829978797172e-01, 1.0732509180989648e-01,  1.0070086522629444e-01,
      9.1016893809887978e-02, 1.0732509180989648e-01,  7.8567420131838608e-02,
      2.8757671678057862e-02, -2.8757671678057862e-02, 1.0070086522629444e-01,
      2.8757671678057862e-02, -6.3730715150116232e-02, -1.1068829978797172e-01,
      9.1016893809887978e-02, -2.8757671678057862e-02, -1.1068829978797172e-01,
      -4.6957584637855494e-02};
  static const double DCTIIb[16] __attribute__((aligned)) = {
      6.3730715150116232e-02,  4.6957584637855494e-02,  2.8757671678057862e-02,
      9.6839714164064644e-03,  -1.0732509180989648e-01, -1.0732509180989648e-01,
      -7.8567420131838608e-02, -2.8757671678057862e-02, 9.6839714164064644e-03,
      9.1016893809887978e-02,  1.0732509180989648e-01,  4.6957584637855494e-02,
      1.0070086522629444e-01,  -9.6839714164064644e-03, -1.0732509180989648e-01,
      -6.3730715150116232e-02};
  static const double DCTIIc[16] __attribute__((aligned)) = {
      6.3730715150116232e-02,  -1.0732509180989648e-01, 9.6839714164064644e-03,
      1.0070086522629444e-01,  4.6957584637855494e-02,  -1.0732509180989648e-01,
      9.1016893809887978e-02,  -9.6839714164064644e-03, 2.8757671678057862e-02,
      -7.8567420131838608e-02, 1.0732509180989648e-01,  -1.0732509180989648e-01,
      9.6839714164064644e-03,  -2.8757671678057862e-02, 4.6957584637855494e-02,
      -6.3730715150116232e-02};
  static const double DCTIId[16] __attribute__((aligned)) = {
      -4.6957584637855494e-02, 1.1068829978797172e-01, -2.8757671678057862e-02,
      -9.1016893809887978e-02, 1.1068829978797172e-01, -6.3730715150116232e-02,
      -2.8757671678057862e-02, 1.0070086522629444e-01, -2.8757671678057862e-02,
      -2.8757671678057862e-02, 7.8567420131838608e-02, -1.0732509180989648e-01,
      -9.1016893809887978e-02, 1.0070086522629444e-01, -1.0732509180989648e-01,
      1.1068829978797172e-01};

  for (size_t i = 0; i < 9; i++) {
    z[i] = input[i * st] + input[(17 - i) * st];
  }
  for (size_t i = 0; i < 4; i++) {
    zpm[i] = z[i] + z[8 - i];
    zpm[5 + i] = z[i] - z[8 - i];
  }
  zpm[4] = z[4];

  for (size_t i = 0; i < 5; i++) {
    const double *t0 = &DCTIIHH0[i * 5];
    out[4 * i] = t0[0] * zpm[0] + t0[1] * zpm[1] + t0[2] * zpm[2] +
                 t0[3] * zpm[3] + t0[4] * zpm[4];
  }
  for (size_t i = 0; i < 4; i++) {
    const double *t1 = &DCTIIHH1[i * 4];
    out[2 + 4 * i] =
        t1[0] * zpm[5] + t1[1] * zpm[6] + t1[2] * zpm[7] + t1[3] * zpm[8];
  }

  for (size_t i = 0; i < 9; i++) {
    z[i] -= 2 * input[(17 - i) * st];
  }

  double z0 = DCTIIa[5] * z[4];
  out[1] = z0;
  out[3] = -z0;
  out[5] = -z0;
  out[7] = z0;
  out[9] = (z[0] - z[1] - z[2] + z[3] + z[4] - z[5] - z[6] + z[7] + z[8]) *
           DCTIIa[5];
  out[11] = -z0;
  out[13] = -z0;
  out[15] = z0;
  out[17] = z0;

  for (size_t i = 0; i < 4; i++) {
    double s0 = 0;
    double s1 = 0;
    for (size_t j = 0; j < 4; j++) {
      s0 += DCTIIa[i * 4 + j] * z[j] + DCTIIb[i * 4 + j] * z[j + 5];
      s1 += DCTIIc[i * 4 + j] * z[j] + DCTIId[i * 4 + j] * z[j + 5];
    }
    out[1 + 2 * i] += s0;
    out[1 + 2 * (i + 5)] += s1;
  }

  for (size_t i = 0; i < 18; i++) {
    output[i * st] = out[i];
  }
  output[0] *= 0.5;
}

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
  const double sqrt_pi = 1.77245385090551602729816e0;
  const size_t N = dplan->N;
  for (size_t i = 0; i < N; i++)
    b[i] = 0.0;
  flops += N;
  if (direction == L2C) {
    const double *a = dplan->a;
    for (size_t n = 0; n < N; n = n + 2) {
      const double *ap = &a[n / 2];
      const double *cp = &u[n];
      const double a0 = ap[0] * M_2_PI;
      for (size_t i = 0; i < N - n; i++) {
        b[i * strides] += a0 * ap[i] * cp[i];
      }
      flops += 3 * (N - n);
    }
    b[0] /= 2;

    flops += N;
  } else {
    double *vn = (double *)fftw_malloc(N * sizeof(double));
    const double *an = dplan->an;
    const double *dn = dplan->dn;
    vn[0] = u[0];
    for (size_t i = 1; i < N; i++)
      vn[i] = u[i * strides] * i;

    for (size_t n = 0; n < N; n++)
      b[n * strides] = sqrt_pi * vn[n] * an[n];

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

///////////
/*
  for (size_t block = 0; block < nL - 1; block++) {
    size_t i0 = block * h;
    for (size_t n = 0; n < h; n = n + 2) {
      const size_t n1 = n / 2;
      const double *ap = &a[i0 + n1];
      const double a0 = a[n1];
      if (strides == 1) {
        double *vp = &output_array[i0];
        const double *up = &input_array[i0 + n];
        for (size_t i = 0; i < h; i++) {
          (*vp++) += a0 * (*ap++) * (*up++);
        }
      } else {
        double *vp = &output_array[i0 * strides];
        const double *up = &input_array[(i0 + n) * strides];
        for (size_t i = 0; i < h; i++) {
          (*vp) += a0 * (*ap++) * (*up);
          vp += strides;
          up += strides;
        }
      }
      flops += h * 3;
    }
  }

  // Last two blocks
  size_t i0 = (nL - 1) * h;
  for (size_t n = 0; n < s; n++) {
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

  for (size_t block = 0; block < nL; block++) {
    size_t i0 = block * h;
    size_t j0 = h + i0;
    const long Nm = lmin(N - j0, h);
    for (size_t n = 0; n < h; n = n + 2) {
      if ((long)(Nm - n) < 0)
        break;
      const size_t n1 = (n + h) / 2;
      const double *ap = &a[i0 + n1];
      const double a0 = a[n1];
      if (strides == 1) {
        double *vp = &output_array[i0];
        const double *up = &input_array[j0 + n];
        for (size_t i = 0; i < (size_t)(Nm - n); i++) {
          (*vp++) += a0 * (*ap++) * (*up++);
        }
      } else {
        double *vp = &output_array[i0 * strides];
        const double *up = &input_array[(j0 + n) * strides];
        for (size_t i = 0; i < (size_t)(Nm - n); i++) {
          (*vp) += a0 * (*ap++) * (*up);
          vp += strides;
          up += strides;
        }
      }
      flops += 3 * (Nm - n);
    }
  }
*/
///////

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
  const double sqrt_pi = 1.77245385090551602729816e0;
  size_t s = fmmplan->s;
  size_t N = fmmplan->N;
  const double *an = fmmplan->dplan->an;
  const double *dn = fmmplan->dplan->dn;
  size_t h = 2 * s;
  size_t nL = N / h;
  size_t flops = 0;
  double *op = &output_array[0];
  const double *ia = &input[0];
  const double *ap = &an[0];
  if (strides == 1) {
    for (size_t i = 0; i < N; i++)
      (*op++) += sqrt_pi * (*ia++) * (*ap++);
  } else {
    for (size_t i = 0; i < N; i++) {
      (*op) += sqrt_pi * (*ia++) * (*ap++);
      op += strides;
    }
  }
  flops += N * 3;

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
          (*vp++) -= d0 * (*ia++) * (*ap++);
        }
      } else {
        double *vp = &output_array[i0 * strides];
        const double *ia = &input[i0 + n];
        for (i = 0; i < lmin(h, 2 * h - n); i++) {
          (*vp) -= d0 * (*ia++) * (*ap++);
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

  // Multiply result by (x+1/2)
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
      for (i = 0; i < m - 2; i = i + 2) {
        zp[i] = x[i] + xp[i];
        zm[i] = x[i] - xp[i];
        zp[i + 1] = x[i + 1] - xp[i + 1];
        zm[i + 1] = x[i + 1] + xp[i + 1];
      }
      zp[i] = x[i] + xp[i];
      zm[i] = x[i] - xp[i];

      const double *a0 = &A[0];
      for (i = 0; i < m - 2; i = i + 2) {
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
        so += (*ap++) * x[j + 1];
      }
      if ((i + m) % 2 == 1)
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

double norm(const double *a, const size_t N) {
  double x = 0;
  for (size_t i = 0; i < N; i++)
    x += a[i] * a[i];
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
    if (plan->B != NULL) {
      fftw_free(plan->B);
      plan->B = NULL;
    }
    if (plan->BT != NULL) {
      fftw_free(plan->BT);
      plan->BT = NULL;
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
  double *a = (double *)fftw_malloc(N * sizeof(double));
  size_t N0 = 256 < N ? 256 : N;
  for (size_t i = 0; i < N0; i++) // Integer, precomputed values
    a[i] = _LambdaI(i);
  size_t DN = 8;
  //size_t maxulp = 0;
  //double maxabs = 0.0;
  //double maxrel = 0.0;
  for (size_t i = N0; i < N; i += DN) {
    a[i] = _LambdaI(i);
    //printf("i %d %d\n", i, maxulp);
    size_t J0 = (size_t)i + DN < N ? (size_t)i + DN : N;
    for (size_t j = i + 1; j < J0; j++) {
      a[j] = a[j - 1] * (1 - 0.5 / j);
      // a[j] = ((j << 1) -1) * a[j - 1] / (j << 1) ;
      //a[j] = a[j - 1] - 0.5 * a[j - 1] / j;

      //double s0 = _LambdaI(j);
      //int ulp = (int)((a[j] - s0) / (nexttoward(a[j], 10) - a[j]));
      //maxulp = (abs(ulp) > maxulp) ? abs(ulp) : maxulp;
      //maxabs = (fabs(a[j] - s0) > maxabs) ? (fabs(a[j] - s0)) : maxabs;
      //maxrel =
      //    (fabs((a[j] - s0) / s0) > maxrel) ? (fabs((a[j] - s0) / s0)) : maxrel;
      //printf("%d  %2.16e  %2.16e %2.16e %2.16e %d %d \n", j, a[j], s0,
      //        fabs(a[j] - s0), fabs((a[j] - s0) / s0), ulp, maxulp);
    }
  }
  //printf("%2.16e %2.16e %d  \n", maxrel, maxabs, maxulp);

  dplan->a = a;
  if (direction == C2L | direction == BOTH) {
    double *dn = (double *)fftw_malloc((N + 1) / 2 * sizeof(double));
    double *an = (double *)fftw_malloc(N * sizeof(double));
    dn[0] = 0;
    an[0] = M_2_SQRTPI; //0.88622692545275801364908374167e0;
    for (size_t i = 1; i < N; i++) {
      an[i] = 1 /(a[i] * (2 * i * i + i)); // Using _Lambda(i-0.5) = 1/(i*_Lambda(i))
    }
    for (size_t i = 1; i < (N + 1) / 2; i++)
      dn[i] = a[i - 1] / (2 * i);

    dplan->an = an;
    dplan->dn = dn;
  }
  dplan->direction = direction;
  dplan->N = N;
  return dplan;
}

fmm_plan *create_fmm(const size_t N, const size_t maxs, const size_t M,
                     const size_t direction, size_t lagrange, const size_t v) {
  fmm_plan *fmmplan = (fmm_plan *)malloc(sizeof(fmm_plan));
  fftw_plan plan1d, plan;
  uint64_t t1 = tic;
  size_t Nn;
  size_t s;
  size_t ij[2];
  size_t directions[2];
  size_t num_directions = 2;
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
  fmmplan->B = NULL;
  fmmplan->BT = NULL;
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
    printf("Lagrange %lu\n", lagrange);
  }

  double *fun = (double *)fftw_malloc(M * M * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(M * M * sizeof(double));
  bool use_FFTW = (M != 18 && lagrange == 0);

  if (use_FFTW) {
    if (v > 1)
      printf("using FFTW for planning\n");
    plan1d = fftw_plan_r2r_1d(M, fun, fun_hat, FFTW_REDFT10, FFTW_PATIENT);
    plan = fftw_plan_r2r_2d(M, M, fun, fun_hat, FFTW_REDFT10, FFTW_REDFT10,
                            FFTW_PATIENT);
  }

  const size_t MM = M * M;
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
  double *xjh = (double *)fftw_malloc(M * sizeof(double));
  for (size_t i = 0; i < M; i++) {
    xj[i] = cos((i + 0.5) * M_PI / M);
  }

  double *fx0 = (double *)fftw_malloc(2 * MM * sizeof(double));
  double *fx1 = (double *)fftw_malloc(MM * sizeof(double));
  double *lx1 = (double *)fftw_malloc(MM * sizeof(double));

  // uint64_t t00 = 0;
  // uint64_t t11 = 0;

  size_t kk = 0;
  for (size_t level = 0; level < L; level++) {
    size_t h = s * get_h(level, L);
    for (size_t k = 0; k < M; k++)
      xjh[k] = xj[k] * h;
    for (size_t block = 0; block < get_number_of_blocks(level); block++) {
      get_ij(ij, level, block, s, L);
      for (size_t q = 0; q < 2; q++) {
        size_t y0 = 2 * (ij[1] + q * h) + h;
        for (size_t p = 0; p < q + 1; p++) {
          size_t x0 = 2 * (ij[0] + p * h) + h;
          for (size_t di = 0; di < num_directions; di++) {
            const size_t dir = directions[di];
            // uint64_t r0 = tic;
            //  ff is input to the DCT
            double *ff = (lagrange == 0) ? &fun[0] : &A[dir][kk * MM];
            double *f00 = &fx0[q * MM];
            double *fpq = &fx0[(q - p) * MM];

            for (size_t i = 0; i < M; i++) {
              double x = x0 + xjh[i];
              for (size_t j = 0; j < M; j++) {
                double y = y0 + xjh[j];
                size_t ix = i * M + j;
                size_t xi = j * M + i;
                if (di == 0 && block == 0 && p == 0) {
                  if (j < M - i) { // Persymmetric Lambda((y-x)/2)
                    double m0 = _Lambda((y - x) / 2);
                    //double m0 = _LambdaE((y - x) / 2);
                    f00[ix] = m0;
                    f00[MM - xi - 1] = m0;
                  }
                }
                if (di == 0) {
                  if (j >= i) { // Symmetric Lambda((x+y)/2)
                    double m1 = _Lambda((x + y) / 2);
                    //double m1 = _LambdaE((x + y) / 2);
                    fx1[ix] = m1;
                    fx1[xi] = m1;
                  }
                }
                if (dir == L2C) {
                  (*ff++) = fpq[ix] * fx1[ix];
                } else {
                  //(*ff++) = -2 * (fpq[ix] / fx1[ix]) /
                  //          ((y - x - 1) * (x + y) * (x + y + 1));
                  (*ff++) = -2 * (fpq[ix] /
                            (fx1[ix] * (x + y) * (x + y + 1) * (y - x - 1)));
                }
              }
            }

            // uint64_t r1 = tic;
            // t00 += r1 - r0;
            if (lagrange == 0) {
              if (use_FFTW) {
                fftw_execute(plan);
                for (size_t i = 0; i < M; i++) {
                  for (size_t j = 0; j < M; j++) {
                    fun_hat[i * M + j] /= (M * M);
                  }
                }
                for (size_t i = 0; i < M; i++)
                  fun_hat[i] /= 2;
                for (size_t i = 0; i < M; i++)
                  fun_hat[i * M] /= 2;
                memcpy(&A[dir][kk * MM], &fun_hat[0], MM * sizeof(double));
              } else {
                dct2(&fun[0], &A[dir][kk * MM]);
              }
            }
            // t11 += tic - r1;
          }
          kk += 1;
        }
      }
    }
  }
  // printf("Time 0: %2.6e   Time 1: %2.6e\n", t00 / 1.0E9, t11 / 1.0E9);

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

  double *B = NULL;
  double *BT = NULL;
  if (!use_FFTW && lagrange == 0) {
    if (v > 1)
      printf("Using exact binomial matrix\n");
    B = (double *)&BMe[0];
    BT = (double *)&BMeT[0];
  } else {
    if (lagrange == 0) {
      B = (double *)fftw_malloc((MM + M) / 2 * sizeof(double));
      BT = (double *)fftw_malloc((MM + M) / 2 * sizeof(double));
      double *Ba = (double *)fftw_malloc(MM * sizeof(double));
      double *BTa = (double *)fftw_malloc(MM * sizeof(double));
      double *th = &Ba[0];
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
      th = &Ba[0];
      double *tht = &BTa[0];
      for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
          tht[i * M + j] = th[j * M + i];
        }
      }
      /// Move to compact storage
      th = &B[0];
      for (size_t k = 0; k < M; k++) {
        for (size_t j = 0; j < k + 1; j++) {
          (*th++) = Ba[k * M + j];
        }
      }
      tht = &BT[0];
      for (size_t k = 0; k < M; k++) {
        for (size_t j = k; j < M; j++) {
          (*tht++) = BTa[k * M + j];
        }
      }
      fftw_free(Ba);
      fftw_free(BTa);

    } else {
      B = (double *)fftw_malloc(2 * MM * sizeof(double));
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
          B[i * M + j] = s0;
          B[MM + i * M + j] = s1;
          sum0 += s0;
          sum1 += s1;
        }
        for (size_t j = 0; j < M; j++) {
          B[i * M + j] /= sum0;
          B[MM + i * M + j] /= sum1;
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
  fftw_free(xjh);
  if (lagrange == 1) {
    fftw_free(wj);
  }
  fmmplan->B = B;
  fmmplan->BT = BT;
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
  double *B = fmmplan->B;
  double *BT = fmmplan->BT;
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
          matvectri(&B[0], wq, &w1[(b0 * 2 + q0) * M], fmmplan->work, M, false);
          flops += MM; //+2*M;
        } else {
          cblas_dgemv(CblasRowMajor, CblasTrans, M, M, 1, &B[0], M, &wq[0], 1,
                      0, &w1[(b0 * 2 + q0) * M], 1);
          cblas_dgemv(CblasRowMajor, CblasTrans, M, M, 1, &B[MM], M, &wq[M], 1,
                      1, &w1[(b0 * 2 + q0) * M], 1);

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
          matvectri(&BT[0], &c0[block * M], &c1[block * 2 * M], NULL, M, true);
          flops += MM;
        } else {
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &B[0], M,
                      &c0[block * M], 1, 1, &c1[block * 2 * M], 1);
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &B[MM], M,
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
