
#ifndef PHARE_TEST_DIAGNOSTIC_FUNC
#define PHARE_TEST_DIAGNOSTIC_FUNC

namespace PHARE_test
{
double density(double x)
{
    return x * x + 2.;
}

double vx(double x)
{
    (void)x;
    return 1.;
}

double vy(double x)
{
    (void)x;
    return 1.;
}

double vz(double x)
{
    (void)x;
    return 1.;
}

double vthx(double x)
{
    (void)x;
    return 1.;
}

double vthy(double x)
{
    (void)x;
    return 1.;
}

double vthz(double x)
{
    (void)x;
    return 1.;
}

double bx(double x)
{
    return 4. * x;
}

double by(double x)
{
    return 5. * x;
}

double bz(double x)
{
    return 6. * x;
}

double ex(double x)
{
    return x;
}

double ey(double x)
{
    return 2. * x;
}

double ez(double x)
{
    return 3. * x;
}

} // namespace PHARE_test

#endif /*PHARE_TEST_DIAGNOSTIC_FUNC*/