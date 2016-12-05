/* Test program for homework 5

 D.R. Reynolds
 Math 6321 @ SMU
 Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "lmm.hpp"

using namespace std;


// Define classes to compute the ODE RHS function and its Jacobian

//    ODE RHS function class -- instantiates a RHSFunction
class BrusselatorRHS: public RHSFunction {
public:
  double a=1.0;
  double b=3.5;
  double ep=5.0e-6;
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = a - (y[2]+1.0)*y[0] + y[0]*y[0]*y[1];
    f[1] = y[2]*y[0] - y[0]*y[0]*y[1];
    f[2] = (b-y[2])/ep - y[2]*y[0];
    return 0;
  }
};

//    ODE RHS Jacobian function class -- instantiates a RHSJacobian
class BrusselatorJac: public RHSJacobian {
public:
  double ep=5.0e-6;
  int Evaluate(double t, vector<double>& y, Matrix& J) {    // evaluates the RHS Jacobian, J(t,y)
    J(0,0) = -(y[2]+1.0) + 2.0*y[0]*y[1];
    J(0,1) = y[0]*y[0];
    J(0,2) = -y[0];
    J(1,0) = y[2] - 2.0*y[0]*y[1];
    J(1,1) = -y[0]*y[0];
    J(1,2) = y[0];
    J(2,0) = -y[2];
    J(2,1) = 0.0;
    J(2,2) = -1.0/ep - y[0];
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




// main routine
int main() {

  // set up test problem
  BrusselatorRHS rhs;
  BrusselatorJac Jac;
  double h = 1e-2;
  vector<double> y0 = y0_Brusselator();
  double t0 = 0.0;
  double Tf = 0.03;
  double tcur = t0;
  Matrix yref = yref_Brusselator();
  cout << "\nRunning Brusselator problem with h = " << h << endl;


  // create BDF method
  vector<double> AM1_a = {1.0};        // aka Trapezoid; used to generate initial conditions
  vector<double> AM1_b = {0.5, 0.5};
  LMMStepper AM1(rhs, Jac, y0, AM1_a, AM1_b);
  vector<double> BDF2_a = {4.0/3.0, -1.0/3.0};
  vector<double> BDF2_b = {2.0/3.0, 0.0, 0.0};
  LMMStepper BDF2(rhs, Jac, y0, BDF2_a, BDF2_b);

  // set Newton solver parameters for AM and BDF methods
  AM1.newt->SetTolerances(1.e-9, 1.e-11);
  AM1.newt->SetMaxit(20);
  AM1.newt->SetShowIterates(true);
  BDF2.newt->SetTolerances(1.e-9, 1.e-11);
  BDF2.newt->SetMaxit(20);
  BDF2.newt->SetShowIterates(true);


  ////////// BDF-2 //////////
  cout << "\n BDF-2:\n";

  // BDF-2 requires two initial conditions
  int M = 3;
  Matrix y(M,2);
  y[1] = y0;        // copy y0 into the 2nd column of y
  //     use AM1 (Trapezoid) to compute second IC
  Matrix ytmp(M,1);
  ytmp[0] = y0;
  vector<double> tspan = {tcur, tcur + h};
  vector<double> tvals = AM1.Evolve(tspan, h, ytmp);
  y[0] = ytmp[0];   // copy ycur into the 1st column of y
  tcur += h;

  // call the solver for the remainder of the time interval
  tspan = {tcur, Tf};
  tvals = BDF2.Evolve(tspan, h, y);
  tcur = tvals.back();      // last entry in tvals

  // compute the error at tcur, output to screen and accumulate maximum
  vector<double> yerr = y[0] - yref[100];
  cout << "\t error = " << InfNorm(yerr) << endl;

  return 0;
}
