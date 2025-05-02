#include <iostream>
#include "capd/capdlib.h"
#include "capd/map/Map.hpp"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/PoincareMap.hpp"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"
#include "capd/dynset/C1DoubletonSet.hpp"
#include <chrono>

using namespace capd;
using capd::autodiff::Node;
using namespace std;

void
cr3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/) {
    // Try to factorize expression as much as possible.
    // Usually we have to define some intermediate quantities.
    Node mj = params[0] - 1; // relative mass of the Sun
    Node xMu = in[0] + params[0];
    Node xMj = in[0] + mj;
    Node xMuSquare = xMu ^ 2; // square
    Node xMjSquare = xMj ^ 2;
    Node ySquare = in[1] ^ 2;
    Node zSquare = in[2] ^ 2;
    Node yzSquare = ySquare + zSquare;

    // power -1.5, for rigorous computation use ONLY REPRESENTABLE CONSTANTS.
    // If exponent is not representable or it is an interval then it should be a parameter of the map.
    Node factor1 = mj * ((xMuSquare + yzSquare) ^ -1.5);
    Node factor2 = params[0] * ((xMjSquare + yzSquare) ^ -1.5);
    Node factor = factor1 - factor2;
    out[0] = in[3];
    out[1] = in[4];
    out[2] = in[5];

    out[3] = in[0] + xMu * factor1 - xMj * factor2 + 2 * in[4];
    out[4] = in[1] * (1 + factor) - 2 * in[3];
    out[5] = in[2] * factor;
}

void g_(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/) {
    Node x2 = sqr(in[0]);
    out[0] = x2 + sqr(in[1]) - sqr(params[0]);
    out[1] = in[1] - x2;
}

DVector newton(DVector x0, DMap &f) {
    for (int i = 0; i < 10; ++i) {
        DMatrix d = f.derivative(x0);
        x0 = x0 - matrixAlgorithms::gauss(d, f(x0));
    }
    return x0;
}

bool intervalNewton(IVector x0, IMap &f, IVector &delta) {
    IVector X = x0 + delta; //sparametryzaować 0.01
    cout << "X:\n" << X << '\n';
    IVector y = f(x0);
    IMatrix d = f.derivative(X);//to już jest otoczka wypukła

    IVector inverseVector = -matrixAlgorithms::gauss(d, y);
    IVector N = x0 + inverseVector;
    cout << "Inverse vector: \n" << inverseVector << '\n';
    cout << "N:\n" << N << '\n';
    if (subset(N, X)) {
        delta = 1.5 * interval(-1, 1) * inverseVector;
        std::cout << "Przybliżenie metodą Newtona: " << x0 << "znajduje się w przedziale:" << N << "\n";
        return true;
    } else {
        std::cout << x0 << "nie jest zerem funkcji." << "\n";
        return false;
    }
}

void calculateSecant() {
    DMap f("var:x,y,a;fun:x^2+y^2-a^2,y-x^2;");
    DVector v{1, 0, 1};
    DMatrix d = f.derivative(v);
    cout << d << '\n';
    DMatrix derivative_u(d.numberOfRows(), d.numberOfColumns() - 1);
    //można jakoś ładnie przekopiować/rozdzielić te macierze?
    for (int i = 0; i < d.numberOfColumns() - 1; ++i) {
        for (int j = 0; j < d.numberOfRows(); ++j) {
            derivative_u[i][j] = d[i][j];
        }
    }
    cout << derivative_u << '\n';
    DVector derivative_a(d.numberOfRows());
    for (int i = 0; i < d.numberOfRows(); ++i) {
        derivative_a[i] = d[i][d.numberOfColumns()];
    }
    cout << derivative_a << '\n';
    DVector secant = matrixAlgorithms::gauss(derivative_u, derivative_a); //double???
    cout << secant << '\n';
}

void newtonIntervalExample() {
    DMap f("var:x,y;par:a;fun:x^2+y^2-a^2,y-x^2;");
    f.setParameter("a", 1);
    IMap g("var:x,y;par:a;fun:x^2+y^2-a,y-x^2;");
//    IMap g(g_,2,2,1);
    g.setParameter(0, 1 + interval(-1, 1) * 0.001);
    DVector x0{0.791667, 0.625};
    x0 = newton(x0, f);
    //intervalNewton(vectalg::convertObject<IVector>(x0), g);

}

void proofOfExistenceOfZeroCurveForInterval() {
    std::cout.precision(17);
    cout << 1 + interval(-1, 1) * 0.001 << '\n';
    //newtonIntervalExample();
    //calculateSecant();
    DMap f("var:x,y;par:a;fun:x^2+y^2-a^2,y-x^2;");
    IMap g("var:x,y;par:a;fun:x^2+y^2-a^2,y-x^2;");

    interval a(1, 2);
    int subintervals = 1;
    double subintervalLength = (a.rightBound() - a.leftBound()) / subintervals;
    DVector point{0.791667, 0.625};
    IVector delta = 0.03 * interval(-1, 1) * IVector{1, 2};
    double left = a.leftBound();
    int counter = 0;
    while (left < a.rightBound()) {
        double right = left + subintervalLength;
        cout << "interval: " << interval(left, right) << '\n';
        g.setParameter("a", interval(left, right));
        double mid = left + ((right - left) / 2.0);// w jakim punkcie zgadujemy 0? w środku przedziału?
        f.setParameter("a", mid);
        DVector x = point;
        x = newton(x, f);
        cout << "mid point: " << x << "value at f: " << f(x) << '\n';
        if (!intervalNewton(vectalg::convertObject<IVector>(x), g, delta)) {
            cout << "Failed for " << left << "\n";
            subintervalLength *= 0.95;
            if (subintervalLength < 1e-6) {
                break;
            }
        } else {
            counter++;
            left = right;
            point = x;
        }
    }
    cout << "subinterval: " << subintervalLength << "steps: " << counter;
}

LDVector findApproxOrbit(LDVector v, LDPoincareMap &pm, int p = 0, int q = 4) {
    for (int i = 0; i < 15; ++i) {
        LDMatrix D(6, 6);
        LDVector y = pm(v, D);
        D = pm.computeDP(y, D);

        // we want y=0 and dx=0
        LDVector u({y[1], y[3]});
        // input variables are x,dy (0 and 4)
        LDMatrix M({{D[1][p], D[1][q]},
                    {D[3][p], D[3][q]}});
        //cout << "DET= " << (M[0][0] * M[1][1] - M[0][1]*M[1][0]) << '\n';
        cout << "u: " << u << '\n';
        u = matrixAlgorithms::gauss(M, u);
        v[p] -= u[0];
        v[q] -= u[1];
        if (u.euclNorm() < 1e-11) {
            break;
        }
    }
    return v;
}

bool intervalNewton(IVector x0, IPoincareMap &pm, IVector &delta, int p = 0, int q = 4) {
    IVector X = x0 + delta;
//    cout << "X:\n" << X << '\n';
    IMatrix D(6, 6);
    C0Rect2Set s0(x0); //0 oznacza, że liczymy tylko przecięcie z sekcją
    IVector y0 = pm(s0);

    C1Rect2Set s(X); //tu liczymy też pochodne odwzorowania poincarego
    IVector y = pm(s, D);
    D = pm.computeDP(y, D);

    IVector u({y0[1], y0[3]});
    IMatrix M({{D[1][p], D[1][q]},//tu 0 i 4 => p i q
               {D[3][p], D[3][q]}});

    IVector inverseVector = -matrixAlgorithms::gauss(M, u);
    IVector N = IVector{x0[p], x0[q]} + inverseVector;
//    cout << "Inverse vector: \n" << inverseVector << '\n';
//    cout << "N:\n" << N << '\n';
    if (subset(N[0], X[p]) && subset(N[1], X[q])) {
//        delta = 1.5 * interval(-1, 1) * inverseVector;
//        cout<<delta<<'\n';
//        delta[0] = 1.1 * interval(-1, 1) * inverseVector[0];
//        delta[4] = 1.1 * interval(-1, 1) * inverseVector[1];
//        cout<<delta<<'\n';
        std::cout << "Ok :)\n";
        return true;
    } else {
        std::cout << "nie :(\n" << "\n";
        return false;
    }
}

void proofOfExistenceOfOneVerticalLyapunovOrbit() {
    cout.precision(17);
    int dim = 6, noParam = 1;
    LDMap vf(cr3bpVectorField, dim, dim, noParam);
    // set value of parameters mu, which is relative mass of Jupiter
    // 0 is index of parameter, 0.0009537 is its new value
    long double mu = 0.0009537;
    vf.setParameter(0, mu);

    // define solver, section and Poincare map
    LDOdeSolver solver(vf, 8);
    LDCoordinateSection section(6, 2); // nowa sekcja Poincarego: z = 0
    LDPoincareMap pm(solver, section);

    LDVector L1{0.98441021177268586, 0, 0.121, 0, -0.039710378195527894, 0};
    LDVector v = findApproxOrbit(L1, pm);

    IMap ivf(cr3bpVectorField, dim, dim, noParam);
    ivf.setParameter(0, mu);
    IOdeSolver isolver(ivf, 10);
    ICoordinateSection isection(6, 2); // nowa sekcja Poincarego: z = 0
    IPoincareMap ipm(isolver, isection);
    IVector delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    intervalNewton(vectalg::convertObject<IVector>(v), ipm, delta);
}

void investigateCurveWhenZIsParametrized() {
    int dim = 6, noParam = 1;
    LDMap vf(cr3bpVectorField, dim, dim, noParam);
    // set value of parameters mu, which is relative mass of Jupiter
    // 0 is index of parameter, 0.0009537 is its new value
    long double mu = 0.0009537;
    vf.setParameter(0, mu);

    // define solver, section and Poincare map
    LDOdeSolver solver(vf, 20);
    LDCoordinateSection section(6, 2); // nowa sekcja Poincarego: z = 0
    LDPoincareMap pm(solver, section);

    LDVector L1{0.87260055880273458, 0, 0.501, 0, 0.09900205284302278, 0};

    std::ofstream outFile("orbits.txt");
    if (!outFile.is_open()) {
        std::cerr << "Nie można otworzyć pliku do zapisu!" << std::endl;
    }
    outFile.precision(17);

    // Pętla zwiększająca trzecią współrzędną (indeks 2) co 1/100
    for (int i = 0; i < 369; ++i) { //wywala sie na 84
        cout << i << '\n';
        long double newValue = L1[2] + i * 0.000001;// zagescic tutaj, popatrzec na skoki
        L1[2] = newValue;

        LDVector v = findApproxOrbit(L1, pm);

        // Zapisanie tylko niezerowych współrzędnych (1, 3, 5 - indeksy 0, 2, 4)
        outFile << v[0] << " "
                << v[2] << " "
                << v[4] << "\n";
    }
    outFile.close();
    std::cout << "Dane zostały zapisane do pliku orbits.txt" << std::endl;
}

void investigateCurveWhenXIsParametrized() {
    int dim = 6, noParam = 1;
    LDMap vf(cr3bpVectorField, dim, dim, noParam);
    // set value of parameters mu, which is relative mass of Jupiter
    // 0 is index of parameter, 0.0009537 is its new value
    long double mu = 0.0009537;
    vf.setParameter(0, mu);

    // define solver, section and Poincare map
    LDOdeSolver solver(vf, 20);
    LDCoordinateSection section(6, 2); // nowa sekcja Poincarego: z = 0
    LDPoincareMap pm(solver, section);

    LDVector L1{-0.99942165641846712, 0, 0.099934538083951164, 0, 1.984734632799688, 0};
//-0.99942165641846712 0.099934538083951164 1.984734632799688
    std::ofstream outFile("orbits.txt");
    if (!outFile.is_open()) {
        std::cerr << "Nie można otworzyć pliku do zapisu!" << std::endl;
    }
    outFile.precision(17);
    double delta = 0.0001;
    long double oldValue = L1[0];
    try {
        while (L1[2] > 0.001) {
            oldValue = L1[0];
            cout << oldValue << '\n';
            long double newValue = L1[0] - delta;
            L1[0] = newValue;

            LDVector v = findApproxOrbit(L1, pm, 2);
            if((v - L1).euclNorm() > 100 * delta) {
                cout << "v- L1 " << (v - L1) << " norma: " << (v - L1).euclNorm() << '\n';
                delta *= 0.9;
                cout << "delta " << delta <<'\n';
                if(delta < 1e-10) {
                    cout << "delta to small"<<'\n';
                    break;
                }
                cout<< "Reducing delta"<<'\n';
                L1[0] = oldValue;
                continue;
            }
            L1 = v;
            delta *= 1.01;
            // Zapisanie tylko niezerowych współrzędnych (1, 3, 5 - indeksy 0, 2, 4)
            outFile << v[0] << " "
                    << v[2] << " "
                    << v[4] << "\n";
        }
    } catch (...){
        outFile.close();
    }
    std::cout << "Dane zostały zapisane do pliku orbits.txt" << std::endl;
}

void proofOfLyapunovOrbits() {
    auto start = std::chrono::high_resolution_clock::now();
    cout.precision(17);
    int dim = 6, noParam = 1;
    LDMap vf(cr3bpVectorField, dim, dim, noParam);
    // set value of parameters mu, which is relative mass of Jupiter
    // 0 is index of parameter, 0.0009537 is its new value
    long double mu = 0.0009537;
    vf.setParameter(0, mu);

    // define solver, section and Poincare map
    LDOdeSolver solver(vf, 20);
    LDCoordinateSection section(6, 2); // nowa sekcja Poincarego: z = 0
    LDPoincareMap pm(solver, section);

//    LDVector L1{0.93236975242513466, 0, 1.3303999999985158e-06, 0, -6.8042716343228538e-12, 0};
//0.87250055880273457 0.50117224951931282 0.099107235220833931
    LDVector L1{0.87250055880273457, 0, 0.50117224951931282, 0, 0.099107235220833931, 0};

    LDVector v = findApproxOrbit(L1, pm);

    IMap ivf(cr3bpVectorField, dim, dim, noParam);
    ivf.setParameter(0, mu);
    IOdeSolver isolver(ivf, 20);
    ICoordinateSection isection(6, 2); // nowa sekcja Poincarego: z = 0
    IPoincareMap ipm(isolver, isection);
    IVector delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    intervalNewton(vectalg::convertObject<IVector>(v), ipm, delta);

    interval a(0.50117224951931282, 0.625);// w strone 0 dojsc do jakiejs malej liczby reprezentowalnej 1/2^20, w prawą stronę do 0.5
    1./(1<<20);//2^20
    int subintervals = 100;
    double subintervalLength = (a.rightBound() - a.leftBound()) / subintervals;
    double left = a.leftBound();
    int counter = 0;
    while (left < a.rightBound()) {
        double right = left + subintervalLength;
        cout << "interval: " << interval(left, right) << '\n';
        double mid = left + ((right - left) / 2.0);
        L1[2] = mid;
        cout<< "L1 candidate: " << L1 << '\n';
        try {
            LDVector x = findApproxOrbit(L1, pm);
            //cout << "mid point: " << x << "value at f: " << f(x) << '\n'
            IVector xInt = vectalg::convertObject<IVector>(x);
            xInt[2] = interval(left, right);
            cout << "Approx orbit: " << xInt << '\n';
            if (!intervalNewton(xInt, ipm, delta)) {
                cout << "Failed for " << left << ". Reducing from " << subintervalLength;
                subintervalLength *= 0.95;
                cout << " to " << subintervalLength << "\n";
                if (subintervalLength < 1e-6) {
                    break;
                }
            } else {
                counter++;
                left = right;
                L1 = x;
                subintervalLength *= 1.01;
            }
        } catch (const std::exception& e) {
            cout << "EXCEPTION for " << left << ". Reducing from " << subintervalLength;
            subintervalLength *= 0.95;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - start;
    cout << "subinterval: " << subintervalLength << "steps: " << counter << " in time: "<<duration.count();
    //subinterval: 1.4989025404881595e-05steps: 506 in time: 3.4105491857000003
    //subinterval: 6.3771611699438274e-06steps: 74076 in time: 475.25548517700003
    //subinterval: 8.8306463463031598e-06steps: 42479 in time: 239.50502835450001 after adding subintervalLength *= 1.01; when success
    //sprobowac pociagnac do 0.625 i sprobowac się przełączyć na x
}

void proofOfLyapunovOrbitsWhenXParametrized() {
    auto start = std::chrono::high_resolution_clock::now();
    cout.precision(17);
    int dim = 6, noParam = 1;
    LDMap vf(cr3bpVectorField, dim, dim, noParam);
    // set value of parameters mu, which is relative mass of Jupiter
    // 0 is index of parameter, 0.0009537 is its new value
    long double mu = 0.0009537;
    vf.setParameter(0, mu);

    // define solver, section and Poincare map
    LDOdeSolver solver(vf, 20);
    LDCoordinateSection section(6, 2); // nowa sekcja Poincarego: z = 0
    LDPoincareMap pm(solver, section);

//    LDVector L1{0.93236975242513466, 0, 1.3303999999985158e-06, 0, -6.8042716343228538e-12, 0};
//0.87250055880273457 0.50117224951931282 0.099107235220833931
    LDVector L1{0.78877760510969482, 0, 0.62372006102301958, 0, 0.18600171965652138, 0};

    LDVector v = findApproxOrbit(L1, pm);

    IMap ivf(cr3bpVectorField, dim, dim, noParam);
    ivf.setParameter(0, mu);
    IOdeSolver isolver(ivf, 20);
    ICoordinateSection isection(6, 2); // nowa sekcja Poincarego: z = 0
    IPoincareMap ipm(isolver, isection);
    IVector delta = IVector{0, 0, 1, 0, 1, 0} * interval(-1, 1) * 1e-5;
    intervalNewton(vectalg::convertObject<IVector>(v), ipm, delta, 2);

    interval a(-1.0045116127308072, 0.78877760510969482);// w strone 0 dojsc do jakiejs malej liczby reprezentowalnej 1/2^20, w prawą stronę do 0.5
    1./(1<<20);//2^20
    int subintervals = 10000;
    double subintervalLength = (a.rightBound() - a.leftBound()) / subintervals;
    double right = a.rightBound();
    double left = a.rightBound();
    int counter = 0;
    while (left > a.leftBound()) {
        left = right - subintervalLength;
        cout << "interval: " << interval(left, right) << '\n';
        double mid = left + ((right - left) / 2.0);
        L1[0] = mid;
        cout<< "L1 candidate: " << L1 << '\n';
        try {
            LDVector x = findApproxOrbit(L1, pm, 2);
            //cout << "mid point: " << x << "value at f: " << f(x) << '\n'
            IVector xInt = vectalg::convertObject<IVector>(x);
            xInt[0] = interval(left, right);
            cout << "Approx orbit: " << xInt << '\n';
            if (!intervalNewton(xInt, ipm, delta, 2)) {
                cout << "Failed for " << left << ". Reducing from " << subintervalLength;
                subintervalLength *= 0.95;
                cout << " to " << subintervalLength << "\n";
                if (subintervalLength < 1e-6) {
                    break;
                }
            } else {
                counter++;
                right = left;
                L1 = x;
                subintervalLength *= 1.01;
            }
        } catch (const std::exception& e) {
            cout << "EXCEPTION for " << left << ". Reducing from " << subintervalLength;
            subintervalLength *= 0.95;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - start;
    cout << "subinterval: " << subintervalLength << "steps: " << counter << " in time: "<<duration.count();
    //subinterval: 1.4989025404881595e-05steps: 506 in time: 3.4105491857000003
    //subinterval: 6.3771611699438274e-06steps: 74076 in time: 475.25548517700003
    //subinterval: 8.8306463463031598e-06steps: 42479 in time: 239.50502835450001 after adding subintervalLength *= 1.01; when success
    //sprobowac pociagnac do 0.625 i sprobowac się przełączyć na x
}

int main() {
    proofOfLyapunovOrbitsWhenXParametrized();
    return 0;
}