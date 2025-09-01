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

void cr3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/) {
    Node mj = params[0] - 1; // relative mass of the Jupiter
    Node xMu = in[0] + params[0];
    Node xMj = in[0] + mj;
    Node xMuSquare = xMu ^ 2;
    Node xMjSquare = xMj ^ 2;
    Node ySquare = in[1] ^ 2;
    Node zSquare = in[2] ^ 2;
    Node yzSquare = ySquare + zSquare;

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
    return make_pair(v, D);
}

auto intervalNewton(IVector x0, IPoincareMap &pm, IVector &delta, const LDMatrix &m, int p = 0, int q = 4) {
    int r = (p == 0) ? 2 : 0;//zakłądamy q stałe = 4, r to ten który parametryzujemy
    LDMatrix DX({{m[1][p], m[1][q]},
                 {m[3][p], m[3][q]}});
    LDVector DZ{m[1][r], m[3][r]};
    LDVector A = matrixAlgorithms::gauss(DX, DZ);//to jest wektor styczny

    IVector x00(6), r0(6);
    split(x0, x00, r0);
    C0Rect2Set s00(x00); //0 oznacza, że liczymy tylko przecięcie z sekcją // g(z0,w0) - złożenie f z s, które w środkowym punkcie daje dokładnie to samo co f
    IVector w0 = pm(s00);
    IVector X = x0 + delta;

    //teraz potrzebujemy policzyć pochodną, na zbiorze Z, w0
    IMatrix DZw0(6, 6);
    IMatrix T = IMatrix::Identity(6);
    T[p][r] = -A[0];
    T[q][r] = -A[1];
    //w tym celu obliczamy C1... (liczy też pochodne odwzorowania) na punkcie x0
    C1Rect2Set s0g = C1Rect2Set(C1Rect2Set::C0BaseSet(x0),
                                C1Rect2Set::C1BaseSet(T)); //tu liczymy też pochodne odwzorowania poincarego
    IVector yg = pm(s0g, DZw0);
    DZw0 = pm.computeDP(yg, DZw0);

    IVector v = w0 + DZw0.column(r) * r0;
    IVector u({v[1], v[3]});
    // g(z0, w0) + [Dzg(Z,w0)]I * Z

    IMatrix DZx0(6, 6);
    split(X, x00, r0);
    C1Rect2Set s0f = C1Rect2Set(C1Rect2Set::C0BaseSet(x00, T, r0),
                                C1Rect2Set::C1BaseSet(T)); //tu liczymy też pochodne odwzorowania poincarego
    IVector yf = pm(s0f, DZx0);
    DZx0 = pm.computeDP(yf, DZx0);

    IMatrix M({{DZx0[1][p], DZx0[1][q]},//tu 0 i 4 => p i q
               {DZx0[3][p], DZx0[3][q]}});

    IVector inverseVector = -matrixAlgorithms::gauss(M, u);
    IVector N = IVector{x0[p], x0[q]} + inverseVector;
    IVector N6d(6);
    N6d[p] = N[0];
    N6d[q] = N[1];
    if (subset(N[0], X[p]) && subset(N[1], X[q])) {
        return make_pair(true, make_pair(N6d, T));
    } else {
        return make_pair(false, make_pair(N6d, T));
    }
}

auto intervalNewton(IVector x0, IPoincareMap &pm, IVector &delta, int p = 0, int q = 4) {
    IVector X = x0 + delta;
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
    IVector N6d(6);
    N6d[p] = N[0];
    N6d[q] = N[1];
    if (subset(N[0], X[p]) && subset(N[1], X[q])) {
        return make_pair(true, N6d);
    } else {
        return make_pair(false, N6d);
    }
}

bool intervalsConnect(int p, int q, const IVector &Nl, const LDVector &x, const IMatrix &T, IVector &delta) {
    IVector Wi = matrixAlgorithms::gauss(T, Nl - vectalg::convertObject<IVector>(x));
    if (subset(Wi[p], delta[p]) && subset(Wi[q], delta[q])) {
        return true;
    } else {
        return false;
    }
}

void proofOfLyapunovOrbitsWithParams(
        interval searchInterval,
        LDPoincareMap &pm,
        IPoincareMap &ipm,
        IVector &delta,
        LDVector L1,
        int parametrizedParam,
        int p,
        bool forward,
        int q = 4
) {
    cout.precision(17);

    int subintervals = 10000;
    double subintervalLength = (searchInterval.rightBound() - searchInterval.leftBound()) / subintervals;
    double left = forward ? searchInterval.leftBound() : searchInterval.rightBound();
    double right = searchInterval.rightBound();

    double L = searchInterval.leftBound();
    double R = searchInterval.rightBound();

    IVector Nl(6);
    bool success = false;
    bool firstStep = true;
    while ((forward && left < searchInterval.rightBound())
           || (!forward && left > searchInterval.leftBound()) || !success) {
        if (subintervalLength < 1e-7) {
            throw std::runtime_error("Interval to small");
        }
        forward ? right = std::min(left + subintervalLength, R)
                : left = std::max(right - subintervalLength, L);
        cout << "\ninterval: " << interval(left, right) << '\n';
        double mid = left + ((right - left) / 2.0);
        L1[parametrizedParam] = mid;

        try {
            auto [x, m] = findApproxOrbit(L1, pm, p);
            auto xInt = vectalg::convertObject<IVector>(x);
            xInt[parametrizedParam] = interval(left, right);

            auto [isZero, pair] = intervalNewton(xInt, ipm, delta, m, p);
            IVector N = pair.first;
            IMatrix T = pair.second;
            IVector X = xInt + delta;

            if (!isZero) {
                subintervalLength *= 0.95;
                success = false;
            } else {
                //check if one site connects
                if (!firstStep && !intervalsConnect(p, q, Nl, x, T, delta)) {
                    throw std::runtime_error("Curves do not connect");
                }
                LDVector v = x;
                v[parametrizedParam] = forward ? right : left;
                v = findApproxOrbit(v, pm, p).first;
                IVector Nz = intervalNewton(vectalg::convertObject<IVector>(v), ipm, delta, p).second;
                Nz[parametrizedParam] = forward ? right : left;
                //check if other site connects
                if (!firstStep && !intervalsConnect(p, q, Nz, x, T, delta)) {
                    throw std::runtime_error("Curves do not connect");
                }
                Nl = Nz;
                forward ? left = right : right = left;
                L1 = x;
                subintervalLength *= 1.01;
                success = true;
                firstStep = false;
            }
        } catch (const std::exception &e) {
            subintervalLength *= 0.95;
            success = false;
        }
    }

}

bool checkIfCurvesConnect(LDVector &L, LDPoincareMap &pm, IPoincareMap &ipm) {
    IVector deltaForNewtonInPointProof = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-6;
    auto x = findApproxOrbit(L, pm, 0).first;
    auto xInt = vectalg::convertObject<IVector>(x);
    auto [isZero, N6d] = intervalNewton(xInt, ipm, deltaForNewtonInPointProof, 0);
    xInt[0] = N6d[0];
    IVector deltaForIntervalsConnectCheck = IVector{0, 0, 1, 0, 1, 0} * interval(-1, 1) * 1e-6;
    isZero = intervalNewton(xInt, ipm, deltaForIntervalsConnectCheck, 2).first;
    return isZero;
}

int main() {
    cout.precision(17);
    auto start = std::chrono::high_resolution_clock::now();
    int dim = 6, noParam = 1;
    LDMap vf(cr3bpVectorField, dim, dim, noParam);
    // set value of parameters mu, which is relative mass of Jupiter
    // 0 is index of parameter, 0.0009537 is its new value
    long double mu = 9.5388114032796904 * 1e-4;
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
    a = interval(1. / (1 << 20), 0.625);
    delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    L1 = LDVector{0.93236545001286, 0, 1. / (1 << 20), 0, -4.68507970920547681e-12, 0}; //guess point
    proofOfLyapunovOrbitsWithParams(a, pm, ipm, delta, L1, 2, 0, true);
    cout << "SWITCHING\n";

    //check if curve is parametrized by line at z = 0.625
    LDVector L{0.78795250494679603, 0, 0.625, 0, 0.18685004432949281, 0};
    if (!checkIfCurvesConnect(L, pm, ipm)) {
        throw std::runtime_error("Curves do not connect");
    }

    //proof when X is parametrized, switched from parametrizing z at z = 0.625
    a = interval(-0.78720775989249825, 0.78795250494679603);
    delta = IVector{0, 0, 1, 0, 1, 0} * interval(-1, 1) * 1e-5;
    L1 = LDVector{0.78795250494679603, 0, 0.625, 0, 0.18685004432949281, 0}; //guess point
    proofOfLyapunovOrbitsWithParams(a, pm, ipm, delta, L1, 0, 2, false);

    //check if curve is parametrized by line at z = 0.625
    L = LDVector{-0.78616953504293408, 0, 0.625, 0, 1.7710684413357896, 0};
    if (!checkIfCurvesConnect(L, pm, ipm)) {
        throw std::runtime_error("Curves do not connect");
    }
    cout << "SWITCHING\n";

    //proof when Z is parametrized, switched from parametrizing x at z = 0.625
    a = interval(1. / (1 << 19), 0.625);
    delta = IVector{1, 0, 0, 0, 1, 0} * interval(-1, 1) * 1e-5;
    L1 = LDVector{-0.78616953504293408, 0, 0.625, 0, 1.7710684413357896, 0}; //guess point
    proofOfLyapunovOrbitsWithParams(a, pm, ipm, delta, L1, 2, 0, false);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - start;
    cout << "time: " << duration.count() << " minutes\n";
    return 0;
}