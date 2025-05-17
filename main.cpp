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

void calculateTangent() {
    DMap f("var:x,y,a;fun:x^2+y^2-a^2,y-x^2;");
    DVector v{1, 0, 1};
    DMatrix d = f.derivative(v);
    cout << d << '\n';
    DMatrix derivative_u(d.numberOfRows(), d.numberOfColumns() - 1);
    //można jakoś ładnie przekopiować/rozdzielić te macierze?
    for (int i = 0; i < d.numberOfColumns() - 1; ++i) {
            derivative_u[i] = d[i];
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
        if (!intervalNewton(vectalg::convertObject<IVector>(x), g, delta )) {
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

auto findApproxOrbit(LDVector v, LDPoincareMap &pm, int p = 0, int q = 4) {
    LDMatrix D(6, 6);
    for (int i = 0; i < 15; ++i) {
        LDVector y = pm(v, D);
        D = pm.computeDP(y, D);

        // we want y=0 and dx=0
        LDVector u({y[1], y[3]});
        // input variables are x,dy (0 and 4)
        LDMatrix M({{D[1][p], D[1][q]},
                    {D[3][p], D[3][q]}});
        u = matrixAlgorithms::gauss(M, u);
        v[p] -= u[0];
        v[q] -= u[1];
        if (u.euclNorm() < 1e-11) {
            break;
        }
    }
    return make_pair(v,D);
}

auto intervalNewton(IVector x0, IPoincareMap &pm, IVector &delta, const LDMatrix &m, int p = 0, int q = 4) {
    int r = (p == 0) ? 2 : 0;//zakłądamy q stałe = 4, r to ten który parametryzujemy
    LDMatrix DX({{m[1][p], m[1][q]},
                 {m[3][p], m[3][q]}});
    LDVector DZ{m[1][r],m[3][r]};
    LDVector A = matrixAlgorithms::gauss(DX,DZ);//to jest 2x2

    IVector x00(6), r0(6);
    split(x0,x00,r0);
    C0Rect2Set s00(x00); //0 oznacza, że liczymy tylko przecięcie z sekcją // g(z0,w0) - złożenie f z s, które w środkowym punkcie daje dokładnie to samo co f
    IVector w0 = pm(s00);
    IVector X = x0 + delta;

    //teraz potrzebujemy policzyć pochodną, na zbiorze Z, w0
    IMatrix DZw0(6, 6);
    IMatrix T = IMatrix::Identity(6);
//    T[r][r] = 1;
    T[p][r] = -A[0];
    T[q][r] = -A[1];
    //w tym celu obliczamy C1... (liczy też pochodne odwzorowania) na punkcie x0
    C1Rect2Set s0g = C1Rect2Set (C1Rect2Set::C0BaseSet(x0),C1Rect2Set::C1BaseSet(T)); //tu liczymy też pochodne odwzorowania poincarego
    IVector yg = pm(s0g, DZw0);
    DZw0 = pm.computeDP(yg, DZw0);

    IVector v = w0 + DZw0.column(r) * r0;
    IVector u({v[1], v[3]});
    // g(z0, w0) + [Dzg(Z,w0)]I * Z

    IMatrix DZx0(6, 6);
    split(X, x00, r0);
    C1Rect2Set s0f = C1Rect2Set (C1Rect2Set::C0BaseSet(x00, T, r0),C1Rect2Set::C1BaseSet(T)); //tu liczymy też pochodne odwzorowania poincarego
    IVector yf = pm(s0f, DZx0);
    DZx0 = pm.computeDP(yf, DZx0);

    IMatrix M({{DZx0[1][p], DZx0[1][q]},//tu 0 i 4 => p i q
               {DZx0[3][p], DZx0[3][q]}});

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
        return make_pair(true, N);
    } else {
        std::cout << "nie :(\n";
        return make_pair(false, N);
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
    LDVector v = findApproxOrbit(L1, pm).first;

    IMap ivf(cr3bpVectorField, dim, dim, noParam);
    ivf.setParameter(0, mu);
    IOdeSolver isolver(ivf, 10);
    ICoordinateSection isection(6, 2); // nowa sekcja Poincarego: z = 0
    IPoincareMap ipm(isolver, isection);
    IVector delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    intervalNewton(vectalg::convertObject<IVector>(v), ipm, delta, LDMatrix(3,2));
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

        LDVector v = findApproxOrbit(L1, pm).first;

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

            LDVector v = findApproxOrbit(L1, pm, 2).first;
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


bool intervalsConnect(IVector Nl, IVector Nr, IVector N) {
//    cout << "\nNl\n" << Nl << "\nNr\n" << Nr  << "\nN\n" << N << "\n";
//    IVector NlNr = vectalg::intersection(Nl, Nr);
//    return subset(N[0], NlNr[0]) && subset(N[1], NlNr[1]);
    return true;
};

void proofOfLyapunovOrbitsWithParams(
        interval searchInterval,
        LDPoincareMap& pm,
        IPoincareMap& ipm,
        IVector& delta,
        LDVector L1,
        int parametrizedParam,
        int p,
        bool forward
) {
    auto start = std::chrono::high_resolution_clock::now();
    cout.precision(17);

    int subintervals = 10000;
    double subintervalLength = (searchInterval.rightBound() - searchInterval.leftBound()) / subintervals;
    double left = forward ? searchInterval.leftBound() : searchInterval.rightBound();
    double right = searchInterval.rightBound();

    cout << "\ninterval: " << interval(left, right) << '\n';
    int counter = 0;

    IVector Nl;
    while ((forward && left < searchInterval.rightBound())
            || (!forward && left > searchInterval.leftBound())) {
        forward ? right = left + subintervalLength
                : left  = right - subintervalLength;
        cout << "\ninterval: " << interval(left, right) << '\n';
        double mid = left + ((right - left) / 2.0);
        L1[parametrizedParam] = mid;
        cout << "L1 candidate: " << L1 << '\n';

//        try {
            auto [x, m] = findApproxOrbit(L1, pm, p);
            IVector xInt = vectalg::convertObject<IVector>(x);
            xInt[parametrizedParam] = interval(left, right);
            cout << "Approx orbit: " << xInt << '\n';

            auto[isZero, Nr] = intervalNewton(xInt, ipm, delta, m, p);
            if (!isZero) {
                cout << "Failed for " << left << ". Reducing from " << subintervalLength;
                subintervalLength *= 0.95;
                cout << " to " << subintervalLength << "\n";
                if (subintervalLength < 1e-6) {
                    break;
                }
            } else {
                if (counter!= 0){
                    LDVector v = x;
                    v[0] = right;
                    v = findApproxOrbit(x, pm, p).first;
                    if (!intervalsConnect(Nl, Nr,intervalNewton(
                            vectalg::convertObject<IVector>(v), ipm, delta,m, p).second)) {
                        cout<< "\nSomething is wrong! \n";
                        break;
                    }
                }
                Nl = Nr;
                forward ? left = right : right = left;
                counter++;
                L1 = x;
                subintervalLength *= 1.01;
            }
//        } catch (const std::exception& e) {
//            cout << "EXCEPTION for " << left << ". Reducing from " << subintervalLength;
//            subintervalLength *= 0.95;
//            cout << " to " << subintervalLength << "\n";
//        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - start;
    cout << "subinterval: " << subintervalLength << " steps: " << counter << " in time: " << duration.count() << " minutes\n";
}



int main() {
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


    IMap ivf(cr3bpVectorField, dim, dim, noParam);
    ivf.setParameter(0, mu);
    IOdeSolver isolver(ivf, 20);
    ICoordinateSection isection(6, 2); // nowa sekcja Poincarego: z = 0
    IPoincareMap ipm(isolver, isection);

    interval a;
    IVector delta;
    LDVector L1;

    //proof when Z is parametrized from as close to zero as possible to z = 0.625
    a = interval(1./(1<<20), 0.625);
    delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    L1 = LDVector{0.93236975242513466, 0, 1.3303999999985158e-06, 0, -6.8042716343228538e-12, 0}; //guess point
    proofOfLyapunovOrbitsWithParams(a, pm, ipm, delta, L1, 2, 0, true);
    cout << "SWITCHING\n";

    //proof when X is parametrized, switched from parametrizing z at z = 0.625
    a = interval(-0.78720775989249825, 0.7883859185128598);
    delta = IVector{0, 0, 1, 0, 1, 0} * interval(-1, 1) * 1e-5;
    L1 = LDVector{0.7883859185128598,0,0.62479468685447181,0,0.18640444496387932,0}; //guess point
    proofOfLyapunovOrbitsWithParams(a, pm, ipm, delta, L1, 0, 2, false);
    cout << "SWITCHING\n";

    //proof when Z is parametrized, switched from parametrizing x at z = ...
    a = interval(1./(1<<10), 0.625);
    delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    L1 = LDVector{0.93236975242513466, 0, 1.3303999999985158e-06, 0, -6.8042716343228538e-12, 0}; //guess point
    proofOfLyapunovOrbitsWithParams(a, pm, ipm, delta, L1, 2, 0, false);
    cout << "SWITCHING\n";

    return 0;
}