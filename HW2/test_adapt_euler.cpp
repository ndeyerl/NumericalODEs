/* Test program for adaptive Forward Euler solvers.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <iomanip>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "adapt_euler.hpp"

using namespace std;


// ODE RHS function class for Brusselator problem
class BrusselatorRHS: public RHSFunction {
public:
  double a=1.0;
  double b=3.5;
  double ep=5.0e-6;
  long int counter=0;
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = a - (y[2]+1)*y[0] + y[0]*y[0]*y[1];
    f[1] = y[2]*y[0] - y[0]*y[0]*y[1];
    f[2] = (b-y[2])/ep - y[2]*y[0];
    counter++;
    return 0;
  }
};

vector<double> y0_Brusselator() { 
  vector<double> y0 = {1.2, 3.1, 3.0};
  return y0;
}

Matrix yref_Brusselator() { 
  Matrix yref(3,101);
  int icol=0;
  yref(0,icol) = 1.200000000000000;  yref(1,icol) = 3.100000000000000;  yref(2,icol++) = 3.000000000000000;
  yref(0,icol) = 1.205372881734043;  yref(1,icol) = 3.074340492475304;  yref(2,icol++) = 3.499978906105430;
  yref(0,icol) = 1.208380541919817;  yref(1,icol) = 3.050623621219999;  yref(2,icol++) = 3.499978853469750;
  yref(0,icol) = 1.208584162242436;  yref(1,icol) = 3.029546593382476;  yref(2,icol++) = 3.499978849903795;
  yref(0,icol) = 1.205561755335091;  yref(1,icol) = 3.011833156313863;  yref(2,icol++) = 3.499978902792281;
  yref(0,icol) = 1.198933355032332;  yref(1,icol) = 2.998205320606491;  yref(2,icol++) = 3.499979018784588;
  yref(0,icol) = 1.188382895080529;  yref(1,icol) = 2.989356174148067;  yref(2,icol++) = 3.499979203411886;
  yref(0,icol) = 1.173677761000278;  yref(1,icol) = 2.985922935409007;  yref(2,icol++) = 3.499979460744984;
  yref(0,icol) = 1.154684301924015;  yref(1,icol) = 2.988462281502648;  yref(2,icol++) = 3.499979793122868;
  yref(0,icol) = 1.131378220265375;  yref(1,icol) = 2.997429497547925;  yref(2,icol++) = 3.499980200970885;
  yref(0,icol) = 1.103849552902046;  yref(1,icol) = 3.013162266929252;  yref(2,icol++) = 3.499980682713558;
  yref(0,icol) = 1.072302695493592;  yref(1,icol) = 3.035869174414269;  yref(2,icol++) = 3.499981234774162;
  yref(0,icol) = 1.037052393869570;  yref(1,icol) = 3.065622450095705;  yref(2,icol++) = 3.499981851644862;
  yref(0,icol) = 0.998516707461129;  yref(1,icol) = 3.102354301161206;  yref(2,icol++) = 3.499982526009847;
  yref(0,icol) = 0.957207634528986;  yref(1,icol) = 3.145856405996426;  yref(2,icol++) = 3.499983248909377;
  yref(0,icol) = 0.913719497416882;  yref(1,icol) = 3.195782695743749;  yref(2,icol++) = 3.499984009943031;
  yref(0,icol) = 0.868714539238278;  yref(1,icol) = 3.251656218427151;  yref(2,icol++) = 3.499984797521759;
  yref(0,icol) = 0.822904757099260;  yref(1,icol) = 3.312881386253573;  yref(2,icol++) = 3.499985599185783;
  yref(0,icol) = 0.777029052106663;  yref(1,icol) = 3.378762951161413;  yref(2,icol++) = 3.499986402004466;
  yref(0,icol) = 0.731825464088198;  yref(1,icol) = 3.448532429031571;  yref(2,icol++) = 3.499987193062192;
  yref(0,icol) = 0.687999520675264;  yref(1,icol) = 3.521381395200213;  yref(2,icol++) = 3.499987960012254;
  yref(0,icol) = 0.646191242632308;  yref(1,icol) = 3.596499383210490;  yref(2,icol++) = 3.499988691654258;
  yref(0,icol) = 0.606944578322248;  yref(1,icol) = 3.673112567167834;  yref(2,icol++) = 3.499989378469026;
  yref(0,icol) = 0.570683443518664;  yref(1,icol) = 3.750518598587906;  yref(2,icol++) = 3.499990013037908;
  yref(0,icol) = 0.537697820820320;  yref(1,icol) = 3.828113299535583;  yref(2,icol++) = 3.499990590286059;
  yref(0,icol) = 0.508141659665723;  yref(1,icol) = 3.905406360829752;  yref(2,icol++) = 3.499991107519204;
  yref(0,icol) = 0.482042148847386;  yref(1,icol) = 3.982025296637137;  yref(2,icol++) = 3.499991564261386;
  yref(0,icol) = 0.459318001338873;  yref(1,icol) = 4.057708982943413;  yref(2,icol++) = 3.499991961934986;
  yref(0,icol) = 0.439803232931858;  yref(1,icol) = 4.132293567460436;  yref(2,icol++) = 3.499992303444615;
  yref(0,icol) = 0.423272714933789;  yref(1,icol) = 4.205694109394249;  yref(2,icol++) = 3.499992592729932;
  yref(0,icol) = 0.409466371021389;  yref(1,icol) = 4.277885065396233;  yref(2,icol++) = 3.499992834342209;
  yref(0,icol) = 0.398109909532515;  yref(1,icol) = 4.348881980718531;  yref(2,icol++) = 3.499993033081506;
  yref(0,icol) = 0.388931067122224;  yref(1,icol) = 4.418725815810870;  yref(2,icol++) = 3.499993193712406;
  yref(0,icol) = 0.381671234168865;  yref(1,icol) = 4.487470499598452;  yref(2,icol++) = 3.499993320760565;
  yref(0,icol) = 0.376092928416716;  yref(1,icol) = 4.555173685130774;  yref(2,icol++) = 3.499993418381920;
  yref(0,icol) = 0.371983887390266;  yref(1,icol) = 4.621890315187126;  yref(2,icol++) = 3.499993490291069;
  yref(0,icol) = 0.369158629616217;  yref(1,icol) = 4.687668443373507;  yref(2,icol++) = 3.499993539733947;
  yref(0,icol) = 0.367458272864130;  yref(1,icol) = 4.752546735309661;  yref(2,icol++) = 3.499993569491003;
  yref(0,icol) = 0.366749266427726;  yref(1,icol) = 4.816553133099824;  yref(2,icol++) = 3.499993581899386;
  yref(0,icol) = 0.366921544490171;  yref(1,icol) = 4.879704258052172;  yref(2,icol++) = 3.499993578885260;
  yref(0,icol) = 0.367886468544388;  yref(1,icol) = 4.942005222096536;  yref(2,icol++) = 3.499993561999812;
  yref(0,icol) = 0.369574812531673;  yref(1,icol) = 5.003449601828070;  yref(2,icol++) = 3.499993532454510;
  yref(0,icol) = 0.371934958558531;  yref(1,icol) = 5.064019394038494;  yref(2,icol++) = 3.499993491152678;
  yref(0,icol) = 0.374931412191157;  yref(1,icol) = 5.123684816846426;  yref(2,icol++) = 3.499993438715480;
  yref(0,icol) = 0.378543710779627;  yref(1,icol) = 5.182403847136872;  yref(2,icol++) = 3.499993375501027;
  yref(0,icol) = 0.382765782235115;  yref(1,icol) = 5.240121394650042;  yref(2,icol++) = 3.499993301615592;
  yref(0,icol) = 0.387605812510198;  yref(1,icol) = 5.296768006584560;  yref(2,icol++) = 3.499993216915938;
  yref(0,icol) = 0.393086696711358;  yref(1,icol) = 5.352257973375798;  yref(2,icol++) = 3.499993121001414;
  yref(0,icol) = 0.399247182609704;  yref(1,icol) = 5.406486663444499;  yref(2,icol++) = 3.499993013193954;
  yref(0,icol) = 0.406143870730815;  yref(1,icol) = 5.459326845854687;  yref(2,icol++) = 3.499992892503073;
  yref(0,icol) = 0.413854321071706;  yref(1,icol) = 5.510623653553934;  yref(2,icol++) = 3.499992757571499;
  yref(0,icol) = 0.422481648971759;  yref(1,icol) = 5.560187676237246;  yref(2,icol++) = 3.499992606594747;
  yref(0,icol) = 0.432161199802329;  yref(1,icol) = 5.607785416917405;  yref(2,icol++) = 3.499992437204323;
  yref(0,icol) = 0.443070223537180;  yref(1,icol) = 5.653125940675673;  yref(2,icol++) = 3.499992246298412;
  yref(0,icol) = 0.455442014899005;  yref(1,icol) = 5.695841882408425;  yref(2,icol++) = 3.499992029794438;
  yref(0,icol) = 0.469586907429109;  yref(1,icol) = 5.735461869169574;  yref(2,icol++) = 3.499991782261680;
  yref(0,icol) = 0.485924125904035;  yref(1,icol) = 5.771369484739742;  yref(2,icol++) = 3.499991496363865;
  yref(0,icol) = 0.505031440634459;  yref(1,icol) = 5.802740432552236;  yref(2,icol++) = 3.499991161990253;
  yref(0,icol) = 0.527725143662305;  yref(1,icol) = 5.828443033985126;  yref(2,icol++) = 3.499990764856090;
  yref(0,icol) = 0.555193973866414;  yref(1,icol) = 5.846874350188411;  yref(2,icol++) = 3.499990284159017;
  yref(0,icol) = 0.589234020932037;  yref(1,icol) = 5.855677429238920;  yref(2,icol++) = 3.499989688468396;
  yref(0,icol) = 0.632684311892872;  yref(1,icol) = 5.851225536428898;  yref(2,icol++) = 3.499988928102888;
  yref(0,icol) = 0.690290991847444;  yref(1,icol) = 5.827615615502185;  yref(2,icol++) = 3.499987920008010;
  yref(0,icol) = 0.770570771022770;  yref(1,icol) = 5.774533431627683;  yref(2,icol++) = 3.499986515147568;
  yref(0,icol) = 0.890270178208291;  yref(1,icol) = 5.672228255712917;  yref(2,icol++) = 3.499984420471558;
  yref(0,icol) = 1.086521530424970;  yref(1,icol) = 5.478032372645838;  yref(2,icol++) = 3.499980986202047;
  yref(0,icol) = 1.455387502973534;  yref(1,icol) = 5.084232970823317;  yref(2,icol++) = 3.499974531360754;
  yref(0,icol) = 2.281768574038728;  yref(1,icol) = 4.176980785112148;  yref(2,icol++) = 3.499960070597381;
  yref(0,icol) = 4.023047481887977;  yref(1,icol) = 2.226932120832307;  yref(2,icol++) = 3.499929599742349;
  yref(0,icol) = 4.990934682333280;  yref(1,icol) = 0.891878606290853;  yref(2,icol++) = 3.499912660888899;
  yref(0,icol) = 4.756544074067366;  yref(1,icol) = 0.735408257385905;  yref(2,icol++) = 3.499916762128852;
  yref(0,icol) = 4.363371458175100;  yref(1,icol) = 0.772482391048185;  yref(2,icol++) = 3.499923642321721;
  yref(0,icol) = 3.983017980282136;  yref(1,icol) = 0.835738650892119;  yref(2,icol++) = 3.499930298252805;
  yref(0,icol) = 3.631013603947149;  yref(1,icol) = 0.907281260805429;  yref(2,icol++) = 3.499936458120017;
  yref(0,icol) = 3.306787293017075;  yref(1,icol) = 0.984839406131653;  yref(2,icol++) = 3.499942131906930;
  yref(0,icol) = 3.007890879743463;  yref(1,icol) = 1.068202218594194;  yref(2,icol++) = 3.499947362450048;
  yref(0,icol) = 2.731794389889210;  yref(1,icol) = 1.157494391378331;  yref(2,icol++) = 3.499952194018858;
  yref(0,icol) = 2.476180738951080;  yref(1,icol) = 1.252871033493724;  yref(2,icol++) = 3.499956667158256;
  yref(0,icol) = 2.239012220358115;  yref(1,icol) = 1.354425918859787;  yref(2,icol++) = 3.499960817524817;
  yref(0,icol) = 2.018560423562559;  yref(1,icol) = 1.462132055067153;  yref(2,icol++) = 3.499964675363101;
  yref(0,icol) = 1.813434143464351;  yref(1,icol) = 1.575781468215588;  yref(2,icol++) = 3.499968265017131;
  yref(0,icol) = 1.622608432952767;  yref(1,icol) = 1.694920965594061;  yref(2,icol++) = 3.499971604421870;
  yref(0,icol) = 1.445446542129438;  yref(1,icol) = 1.818792358239893;  yref(2,icol++) = 3.499974704719190;
  yref(0,icol) = 1.281698930209825;  yref(1,icol) = 1.946294407080820;  yref(2,icol++) = 3.499977570275059;
  yref(0,icol) = 1.131459628420471;  yref(1,icol) = 2.075989485538215;  yref(2,icol++) = 3.499980199443058;
  yref(0,icol) = 0.995064137018778;  yref(1,icol) = 2.206175908776501;  yref(2,icol++) = 3.499982586351075;
  yref(0,icol) = 0.872928432151910;  yref(1,icol) = 2.335032398726637;  yref(2,icol++) = 3.499984723718588;
  yref(0,icol) = 0.765353284831617;  yref(1,icol) = 2.460815381730033;  yref(2,icol++) = 3.499986606281041;
  yref(0,icol) = 0.672340674823181;  yref(1,icol) = 2.582063602797944;  yref(2,icol++) = 3.499988234002641;
  yref(0,icol) = 0.593474229385347;  yref(1,icol) = 2.697754259758687;  yref(2,icol++) = 3.499989614168765;
  yref(0,icol) = 0.527896102562875;  yref(1,icol) = 2.807369965127945;  yref(2,icol++) = 3.499990761790685;
  yref(0,icol) = 0.474377827466704;  yref(1,icol) = 2.910869178808560;  yref(2,icol++) = 3.499991698365740;
  yref(0,icol) = 0.431452159767283;  yref(1,icol) = 3.008585158126683;  yref(2,icol++) = 3.499992449570112;
  yref(0,icol) = 0.397561495452604;  yref(1,icol) = 3.101093960065049;  yref(2,icol++) = 3.499993042661507;
  yref(0,icol) = 0.371186500509693;  yref(1,icol) = 3.189088144162925;  yref(2,icol++) = 3.499993504228089;
  yref(0,icol) = 0.350936145846172;  yref(1,icol) = 3.273278056430489;  yref(2,icol++) = 3.499993858612817;
  yref(0,icol) = 0.335596749263632;  yref(1,icol) = 3.354327184953578;  yref(2,icol++) = 3.499994127055157;
  yref(0,icol) = 0.324147547028770;  yref(1,icol) = 3.432817840402948;  yref(2,icol++) = 3.499994327418549;
  yref(0,icol) = 0.315753667988019;  yref(1,icol) = 3.509239098111851;  yref(2,icol++) = 3.499994474313320;
  yref(0,icol) = 0.309746537216893;  yref(1,icol) = 3.583988715664158;  yref(2,icol++) = 3.499994579439618;
  yref(0,icol) = 0.305599158422400;  yref(1,icol) = 3.657382437597814;  yref(2,icol++) = 3.499994652019957;
  return yref;
}

// ODE RHS function class for VanderPol oscillator problem
class VanderPolRHS: public RHSFunction {
public:
  double ep=1.0e-4;
  double ep2=1.0;
  long int counter=0;
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = y[1];
    f[1] = (y[1] - y[1]*y[0]*y[0] - ep2*y[0])/ep;
    counter++;
    return 0;
  }
};

vector<double> y0_VanderPol() { 
  vector<double> y0 = {-2.0, -2.3553013976081};
  return y0;
}

//    Convenience function for reference solution
Matrix yref_VanderPol() { 
  Matrix yref(2,151);
  int icol=0;
  yref(0,icol) = -2.000000000000000;  yref(1,icol++) = -2.355301397608100;
  yref(0,icol) = -1.993416117230992;  yref(1,icol++) = 0.670334323610993;
  yref(0,icol) = -1.986693808859355;  yref(1,icol++) = 0.674138980200715;
  yref(0,icol) = -1.979933101365236;  yref(1,icol++) = 0.678014518457502;
  yref(0,icol) = -1.973133274655014;  yref(1,icol++) = 0.681963210638489;
  yref(0,icol) = -1.966293585397327;  yref(1,icol++) = 0.685987432988440;
  yref(0,icol) = -1.959413265952723;  yref(1,icol++) = 0.690089671956446;
  yref(0,icol) = -1.952491523238669;  yref(1,icol++) = 0.694272530873933;
  yref(0,icol) = -1.945527537525229;  yref(1,icol++) = 0.698538737135234;
  yref(0,icol) = -1.938520461156245;  yref(1,icol++) = 0.702891149924613;
  yref(0,icol) = -1.931469417190339;  yref(1,icol++) = 0.707332768541346;
  yref(0,icol) = -1.924373497955442;  yref(1,icol++) = 0.711866741377217;
  yref(0,icol) = -1.917231763510038;  yref(1,icol++) = 0.716496375608137;
  yref(0,icol) = -1.910043240003651;  yref(1,icol++) = 0.721225147668066;
  yref(0,icol) = -1.902806917928317;  yref(1,icol++) = 0.726056714581322;
  yref(0,icol) = -1.895521750252147;  yref(1,icol++) = 0.730994926237149;
  yref(0,icol) = -1.888186650425055;  yref(1,icol++) = 0.736043838700683;
  yref(0,icol) = -1.880800490245693;  yref(1,icol++) = 0.741207728666356;
  yref(0,icol) = -1.873362097577702;  yref(1,icol++) = 0.746491109170351;
  yref(0,icol) = -1.865870253902011;  yref(1,icol++) = 0.751898746694351;
  yref(0,icol) = -1.858323691690602;  yref(1,icol++) = 0.757435679808666;
  yref(0,icol) = -1.850721091585188;  yref(1,icol++) = 0.763107239521795;
  yref(0,icol) = -1.843061079363300;  yref(1,icol++) = 0.768919071523526;
  yref(0,icol) = -1.835342222671745;  yref(1,icol++) = 0.774877160533490;
  yref(0,icol) = -1.827563027504794;  yref(1,icol++) = 0.780987856995126;
  yref(0,icol) = -1.819721934402649;  yref(1,icol++) = 0.787257906388971;
  yref(0,icol) = -1.811817314342629;  yref(1,icol++) = 0.793694481470956;
  yref(0,icol) = -1.803847464291444;  yref(1,icol++) = 0.800305217786290;
  yref(0,icol) = -1.795810602384661;  yref(1,icol++) = 0.807098252873195;
  yref(0,icol) = -1.787704862694160;  yref(1,icol++) = 0.814082269595045;
  yref(0,icol) = -1.779528289538344;  yref(1,icol++) = 0.821266544146926;
  yref(0,icol) = -1.771278831287426;  yref(1,icol++) = 0.828660999331762;
  yref(0,icol) = -1.762954333606286;  yref(1,icol++) = 0.836276263811866;
  yref(0,icol) = -1.754552532072065;  yref(1,icol++) = 0.844123738137720;
  yref(0,icol) = -1.746071044092936;  yref(1,icol++) = 0.852215668492825;
  yref(0,icol) = -1.737507360046958;  yref(1,icol++) = 0.860565229250877;
  yref(0,icol) = -1.728858833546636;  yref(1,icol++) = 0.869186615603345;
  yref(0,icol) = -1.720122670719678;  yref(1,icol++) = 0.878095147756921;
  yref(0,icol) = -1.711295918382278;  yref(1,icol++) = 0.887307388442446;
  yref(0,icol) = -1.702375450964369;  yref(1,icol++) = 0.896841275795827;
  yref(0,icol) = -1.693357956017836;  yref(1,icol++) = 0.906716274041775;
  yref(0,icol) = -1.684239918117569;  yref(1,icol++) = 0.916953544879730;
  yref(0,icol) = -1.675017600930353;  yref(1,icol++) = 0.927576143012184;
  yref(0,icol) = -1.665687027191027;  yref(1,icol++) = 0.938609239948444;
  yref(0,icol) = -1.656243956279820;  yref(1,icol++) = 0.950080381045101;
  yref(0,icol) = -1.646683859038184;  yref(1,icol++) = 0.962019781779863;
  yref(0,icol) = -1.637001889395294;  yref(1,icol++) = 0.974460670543179;
  yref(0,icol) = -1.627192852298222;  yref(1,icol++) = 0.987439686840792;
  yref(0,icol) = -1.617251167339095;  yref(1,icol++) = 1.000997345802252;
  yref(0,icol) = -1.607170827350000;  yref(1,icol++) = 1.015178582471659;
  yref(0,icol) = -1.596945351089029;  yref(1,icol++) = 1.030033392593691;
  yref(0,icol) = -1.586567728948702;  yref(1,icol++) = 1.045617590799478;
  yref(0,icol) = -1.576030360391493;  yref(1,icol++) = 1.061993712483819;
  yref(0,icol) = -1.565324981510698;  yref(1,icol++) = 1.079232092728120;
  yref(0,icol) = -1.554442580744579;  yref(1,icol++) = 1.097412164868123;
  yref(0,icol) = -1.543373300284036;  yref(1,icol++) = 1.116624033598689;
  yref(0,icol) = -1.532106320087035;  yref(1,icol++) = 1.136970393968191;
  yref(0,icol) = -1.520629720595194;  yref(1,icol++) = 1.158568889904332;
  yref(0,icol) = -1.508930319165463;  yref(1,icol++) = 1.181555036444568;
  yref(0,icol) = -1.496993473788863;  yref(1,icol++) = 1.206085872156554;
  yref(0,icol) = -1.484802845723389;  yref(1,icol++) = 1.232344567677901;
  yref(0,icol) = -1.472340110010274;  yref(1,icol++) = 1.260546300999598;
  yref(0,icol) = -1.459584599163903;  yref(1,icol++) = 1.290945832674050;
  yref(0,icol) = -1.446512860148530;  yref(1,icol++) = 1.323847394564699;
  yref(0,icol) = -1.433098097360386;  yref(1,icol++) = 1.359617776333641;
  yref(0,icol) = -1.419309463573325;  yref(1,icol++) = 1.398703908095646;
  yref(0,icol) = -1.405111144819044;  yref(1,icol++) = 1.441656886551825;
  yref(0,icol) = -1.390461160913629;  yref(1,icol++) = 1.489165434353729;
  yref(0,icol) = -1.375309765571507;  yref(1,icol++) = 1.542103505957805;
  yref(0,icol) = -1.359597269534539;  yref(1,icol++) = 1.601599696244003;
  yref(0,icol) = -1.343251010070955;  yref(1,icol++) = 1.669141322543229;
  yref(0,icol) = -1.326181018369715;  yref(1,icol++) = 1.746735690427306;
  yref(0,icol) = -1.308273628289912;  yref(1,icol++) = 1.837169780994378;
  yref(0,icol) = -1.289381688634903;  yref(1,icol++) = 1.944448176883524;
  yref(0,icol) = -1.269308874807151;  yref(1,icol++) = 2.074574335348090;
  yref(0,icol) = -1.247783072641715;  yref(1,icol++) = 2.237045987277913;
  yref(0,icol) = -1.224407806199182;  yref(1,icol++) = 2.447989019107501;
  yref(0,icol) = -1.198564525354610;  yref(1,icol++) = 2.737576011806878;
  yref(0,icol) = -1.169186979466559;  yref(1,icol++) = 3.170925674663635;
  yref(0,icol) = -1.134116109086221;  yref(1,icol++) = 3.925580742929147;
  yref(0,icol) = -1.087361638841265;  yref(1,icol++) = 5.778951465443521;
  yref(0,icol) = -0.970981398215555;  yref(1,icol++) = 40.212565136715803;
  yref(0,icol) = 1.996218411070505;  yref(1,icol++) = -0.668762457503965;
  yref(0,icol) = 1.989511965687082;  yref(1,icol++) = -0.672538098432962;
  yref(0,icol) = 1.982767415911700;  yref(1,icol++) = -0.676383699477365;
  yref(0,icol) = 1.975984051073178;  yref(1,icol++) = -0.680301491114734;
  yref(0,icol) = 1.969161137692592;  yref(1,icol++) = -0.684293805336132;
  yref(0,icol) = 1.962297918438600;  yref(1,icol++) = -0.688363081680995;
  yref(0,icol) = 1.955393611019963;  yref(1,icol++) = -0.692511873717438;
  yref(0,icol) = 1.948447407010836;  yref(1,icol++) = -0.696742856007442;
  yref(0,icol) = 1.941458470603789;  yref(1,icol++) = -0.701058831598998;
  yref(0,icol) = 1.934425937285111;  yref(1,icol++) = -0.705462740094792;
  yref(0,icol) = 1.927348912426322;  yref(1,icol++) = -0.709957666349395;
  yref(0,icol) = 1.920226469785353;  yref(1,icol++) = -0.714546849853993;
  yref(0,icol) = 1.913057649910124;  yref(1,icol++) = -0.719233694874329;
  yref(0,icol) = 1.905841458436702;  yref(1,icol++) = -0.724021781413970;
  yref(0,icol) = 1.898576864273375;  yref(1,icol++) = -0.728914877083473;
  yref(0,icol) = 1.891262797661148;  yref(1,icol++) = -0.733916949965494;
  yref(0,icol) = 1.883898148100126;  yref(1,icol++) = -0.739032182576738;
  yref(0,icol) = 1.876481762130394;  yref(1,icol++) = -0.744264987038068;
  yref(0,icol) = 1.869012440954674;  yref(1,icol++) = -0.749620021578554;
  yref(0,icol) = 1.861488937888615;  yref(1,icol++) = -0.755102208515938;
  yref(0,icol) = 1.853909955623240;  yref(1,icol++) = -0.760716753869909;
  yref(0,icol) = 1.846274143282275;  yref(1,icol++) = -0.766469168788370;
  yref(0,icol) = 1.838580093255651;  yref(1,icol++) = -0.772365292987528;
  yref(0,icol) = 1.830826337787304;  yref(1,icol++) = -0.778411320432854;
  yref(0,icol) = 1.823011345293895;  yref(1,icol++) = -0.784613827520710;
  yref(0,icol) = 1.815133516388280;  yref(1,icol++) = -0.790979804052599;
  yref(0,icol) = 1.807191179577477;  yref(1,icol++) = -0.797516687331954;
  yref(0,icol) = 1.799182586602721;  yref(1,icol++) = -0.804232399771111;
  yref(0,icol) = 1.791105907384621;  yref(1,icol++) = -0.811135390435655;
  yref(0,icol) = 1.782959224530045;  yref(1,icol++) = -0.818234681025507;
  yref(0,icol) = 1.774740527355431;  yref(1,icol++) = -0.825539916865465;
  yref(0,icol) = 1.766447705372280;  yref(1,icol++) = -0.833061423567448;
  yref(0,icol) = 1.758078541174592;  yref(1,icol++) = -0.840810270112019;
  yref(0,icol) = 1.749630702660041;  yref(1,icol++) = -0.848798339243332;
  yref(0,icol) = 1.741101734505603;  yref(1,icol++) = -0.857038406193819;
  yref(0,icol) = 1.732489048809690;  yref(1,icol++) = -0.865544226929630;
  yref(0,icol) = 1.723789914799578;  yref(1,icol++) = -0.874330637310051;
  yref(0,icol) = 1.715001447485285;  yref(1,icol++) = -0.883413664791135;
  yref(0,icol) = 1.706120595124256;  yref(1,icol++) = -0.892810654601808;
  yref(0,icol) = 1.697144125345399;  yref(1,icol++) = -0.902540412652925;
  yref(0,icol) = 1.688068609747565;  yref(1,icol++) = -0.912623367872461;
  yref(0,icol) = 1.678890406764614;  yref(1,icol++) = -0.923081757170716;
  yref(0,icol) = 1.669605642553330;  yref(1,icol++) = -0.933939836858452;
  yref(0,icol) = 1.660210189613432;  yref(1,icol++) = -0.945224125108576;
  yref(0,icol) = 1.650699642809226;  yref(1,icol++) = -0.956963681005692;
  yref(0,icol) = 1.641069292385132;  yref(1,icol++) = -0.969190426887756;
  yref(0,icol) = 1.631314093513504;  yref(1,icol++) = -0.981939522155194;
  yref(0,icol) = 1.621428631801470;  yref(1,icol++) = -0.995249798551626;
  yref(0,icol) = 1.611407084090823;  yref(1,icol++) = -1.009164269234066;
  yref(0,icol) = 1.601243173734263;  yref(1,icol++) = -1.023730726894602;
  yref(0,icol) = 1.590930119368326;  yref(1,icol++) = -1.039002449951348;
  yref(0,icol) = 1.580460575986945;  yref(1,icol++) = -1.055039040676252;
  yref(0,icol) = 1.569826566851540;  yref(1,icol++) = -1.071907425427172;
  yref(0,icol) = 1.559019404431633;  yref(1,icol++) = -1.089683055393576;
  yref(0,icol) = 1.548029598136434;  yref(1,icol++) = -1.108451357174664;
  yref(0,icol) = 1.536846746031937;  yref(1,icol++) = -1.128309497050455;
  yref(0,icol) = 1.525459407009090;  yref(1,icol++) = -1.149368542417970;
  yref(0,icol) = 1.513854948906307;  yref(1,icol++) = -1.171756130587540;
  yref(0,icol) = 1.502019366812657;  yref(1,icol++) = -1.195619791990902;
  yref(0,icol) = 1.489937064066956;  yref(1,icol++) = -1.221131126321937;
  yref(0,icol) = 1.477590586136895;  yref(1,icol++) = -1.248491103038526;
  yref(0,icol) = 1.464960294360120;  yref(1,icol++) = -1.277936862424724;
  yref(0,icol) = 1.452023962045692;  yref(1,icol++) = -1.309750546568834;
  yref(0,icol) = 1.438756269082140;  yref(1,icol++) = -1.344270917431013;
  yref(0,icol) = 1.425128162023024;  yref(1,icol++) = -1.381908864850812;
  yref(0,icol) = 1.411106033107567;  yref(1,icol++) = -1.423168443549174;
  yref(0,icol) = 1.396650651357579;  yref(1,icol++) = -1.468675930196526;
  yref(0,icol) = 1.381715747617306;  yref(1,icol++) = -1.519220782883151;
  yref(0,icol) = 1.366246105895326;  yref(1,icol++) = -1.575814728003287;
  return yref;
};

// ODE RHS function class for simple nonlinear problem
class NonlinearRHS: public RHSFunction {
public:
  long int counter=0;
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = (t+1.0)*exp(-y[0]);
    counter++;
    return 0;
  }
};

vector<double> y0_Nonlinear() { 
  vector<double> y0 = {0.0};
  return y0;
}

//    Convenience function for analytical solution
vector<double> ytrue_Nonlinear(const double t) { 
  vector<double> yt(1);
  yt[0] = log(((t+1.0)*(t+1.0)+1.0)/2.0);
  return yt;
};





// main routine
int main() {

  // reusable tolerance arrays
  vector<double> rtols(3);
  double atol;

  // overall statistics arrays
  vector<double> err_tol_ratio(9);  // 3 problems. 3 tolerances
  vector<int> step_attempts(9);
  int itest=0;


  //----- Test problem 1: brusselator -----//
  rtols[0] = 1.e-1;  rtols[1] = 1.e-3;  rtols[2] = 1.e-5;
  atol = 1.e-14;
  for (int itol=0; itol<3; itol++) {
    BrusselatorRHS rhs;
    vector<double> y0 = y0_Brusselator();
    vector<double> y = y0;
    double t0 = 0.0;
    double Tf = 10.0;
    double tcur = t0;
    int Nt = 100;
    double dtout = (Tf-t0)/Nt;
    AdaptEuler AE(rhs, rtols[itol], atol, y0);
    Matrix yref = yref_Brusselator();
    Matrix ycalc = yref;  ycalc[0] = y0;
    double maxerr = 0.0;
    double maxy = InfNorm(yref);
    cout << "\nRunning Brusselator problem with rtol = " << rtols[itol]
         << " and atol = " << atol << endl;

    //   loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {t0+istep*dtout, t0+(istep+1)*dtout};

      // call the solver for this time interval
      vector<double> tvals = AE.Evolve(tspan, y);
      tcur = tvals.back();  // last entry in tvals

      // compute the errors at tcur, output to screen, and accumulate maxima
      ycalc[istep+1] = y;
      vector<double> ytrue = yref[istep+1];
      vector<double> yerr = y - ytrue;
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);
      if (std::abs(tcur - (t0+(istep+1)*dtout))/std::abs(Tf) > 1e-4)
	cout << "  warning: tcur = " << tcur << " differs from " << t0+(istep+1)*dtout << "\n";
    }

    // compute/store stats
    err_tol_ratio[itest] = maxerr / (rtols[itol]*maxy + atol);
    step_attempts[itest] = rhs.counter;

    // output final results for this tolerance
    cout << "  maxerr = " << maxerr << endl
	 << "  goal = " << rtols[itol]*maxy + atol << endl
	 << "  ratio = " << err_tol_ratio[itest] << endl
	 << "  steps = " << step_attempts[itest] << endl;
    itest++;

    // output results to disk
    Matrix tspan = Linspace(t0, Tf, 1, Nt+1);
    tspan.Write("t_Brusselator.txt");
    yref.Write("ytrue_Brusselator.txt");
    char fname[25];
    sprintf(fname, "ycalc_Brusselator%i.txt", itol);
    ycalc.Write(fname);

  }



  //----- Test problem 2: Van der Pol oscillator -----//
  rtols[0] = 1.e-1;  rtols[1] = 1.e-3;  rtols[2] = 1.e-5;
  atol = 1.e-14;
  for (int itol=0; itol<3; itol++) {
    VanderPolRHS rhs;
    vector<double> y0 = y0_VanderPol();
    vector<double> y = y0;
    double t0 = 0.0;
    double Tf = 1.5;
    double tcur = t0;
    int Nt = 150;
    double dtout = (Tf-t0)/Nt;
    AdaptEuler AE(rhs, rtols[itol], atol, y0);
    Matrix yref = yref_VanderPol();
    Matrix ycalc = yref;  ycalc[0] = y0;
    double maxerr = 0.0;
    double maxy = InfNorm(yref);
    cout << "\nRunning Van der Pol problem with rtol = " << rtols[itol]
         << " and atol = " << atol << endl;

    //   loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {t0+istep*dtout, t0+(istep+1)*dtout};

      // call the solver for this time interval
      vector<double> tvals = AE.Evolve(tspan, y);
      tcur = tvals.back();  // last entry in tvals

      // compute the errors at tcur, output to screen, and accumulate maxima
      ycalc[istep+1] = y;
      vector<double> ytrue = yref[istep+1];
      vector<double> yerr = y - ytrue;
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);
      if (std::abs(tcur - (t0+(istep+1)*dtout))/std::abs(Tf) > 1e-4)
	cout << "  warning: tcur = " << tcur << " differs from " << t0+(istep+1)*dtout << "\n";
    }

    // compute/store stats
    err_tol_ratio[itest] = maxerr / (rtols[itol]*maxy + atol);
    step_attempts[itest] = rhs.counter;

    // output final results for this tolerance
    cout << "  maxerr = " << maxerr << endl
	 << "  goal = " << rtols[itol]*maxy + atol << endl
	 << "  ratio = " << err_tol_ratio[itest] << endl
	 << "  steps = " << step_attempts[itest] << endl;
    itest++;

    // output results to disk
    Matrix tspan = Linspace(t0, Tf, 1, Nt+1);
    tspan.Write("t_VanderPol.txt");
    yref.Write("ytrue_VanderPol.txt");
    char fname[25];
    sprintf(fname, "ycalc_VanderPol%i.txt", itol);
    ycalc.Write(fname);

  }



  //----- Test problem 3: simple nonlinear problem -----//
  rtols[0] = 1.e-3;  rtols[1] = 1.e-5;  rtols[2] = 1.e-7;
  atol = 1.e-14;
  for (int itol=0; itol<3; itol++) {
    NonlinearRHS rhs;
    vector<double> y0 = y0_Nonlinear();
    vector<double> y = y0;
    double t0 = 0.0;
    double Tf = 10.0;
    double tcur = t0;
    int Nt = 100;
    double dtout = (Tf-t0)/Nt;
    AdaptEuler AE(rhs, rtols[itol], atol, y0);
    Matrix yref(1,Nt+1);  yref[0] = y0;
    Matrix ycalc = yref;  ycalc[0] = y0;
    double maxerr = 0.0;
    double maxy = 0.0;
    cout << "\nRunning simple nonlinear problem with rtol = " << rtols[itol]
         << " and atol = " << atol << endl;

    //   loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {t0+istep*dtout, t0+(istep+1)*dtout};

      // call the solver for this time interval
      vector<double> tvals = AE.Evolve(tspan, y);
      tcur = tvals.back();  // last entry in tvals

      // compute the errors at tcur, output to screen, and accumulate maxima
      vector<double> ytrue = ytrue_Nonlinear(t0+(istep+1)*dtout);
      vector<double> yerr = y - ytrue;
      double err = InfNorm(yerr);
      yref[istep+1] = ytrue;
      ycalc[istep+1] = y;
      maxerr = std::max(maxerr, err);
      maxy = std::max(maxy, InfNorm(ytrue));
      if (std::abs(tcur - (t0+(istep+1)*dtout))/std::abs(Tf) > 1e-4)
	cout << "  warning: tcur = " << tcur << " differs from " << t0+(istep+1)*dtout << "\n";
    }

    // compute/store stats
    err_tol_ratio[itest] = maxerr / (rtols[itol]*maxy + atol);
    step_attempts[itest] = rhs.counter;

    // output final results for this tolerance
    cout << "  maxerr = " << maxerr << endl
	 << "  goal = " << rtols[itol]*maxy + atol << endl
	 << "  ratio = " << err_tol_ratio[itest] << endl
	 << "  steps = " << step_attempts[itest] << endl;
    itest++;

    // output results to disk
    Matrix tspan = Linspace(t0, Tf, 1, Nt+1);
    tspan.Write("t_Nonlinear.txt");
    yref.Write("ytrue_Nonlinear.txt");
    char fname[25];
    sprintf(fname, "ycalc_Nonlinear%i.txt", itol);
    ycalc.Write(fname);

  }

  // output overall statistics
  double err_tol_overall = 1.0;
  double steps_overall = 1.0;
  for (int i=0; i<itest; i++) {
    err_tol_overall *= pow(err_tol_ratio[i], 1.0/9.0);
    steps_overall *= pow(1.0*step_attempts[i], 1.0/9.0);
  }
  cout << "\n\nFinal Scores:" << "\n"
       << "  error ratio = " << err_tol_overall << "\n"
       << "  work = " << steps_overall << "\n\n";

  return 0;
}
