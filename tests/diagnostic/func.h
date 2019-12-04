
#ifndef PHARE_TEST_DIAGNOSTIC_FUNC
#define PHARE_TEST_DIAGNOSTIC_FUNC

namespace PHARE_test
{
double density(double x)
{
    return x * x + 2.;
}

double vx(double /*x*/)
{
    return 1.;
}

double vy(double /*x*/)
{
    return 1.;
}

double vz(double /*x*/)
{
    return 1.;
}

double vthx(double /*x*/)
{
    return 1.;
}

double vthy(double /*x*/)
{
    return 1.;
}

double vthz(double /*x*/)
{
    return 1.;
}

double bx(double /*x*/)
{
    return 4.;
}

double by(double /*x*/)
{
    return 5.;
}

double bz(double /*x*/)
{
    return 6.;
}

double ex(double /*x*/)
{
    return 1;
}

double ey(double /*x*/)
{
    return 2.;
}

double ez(double /*x*/)
{
    return 3.;
}

} // namespace PHARE_test

#endif /*PHARE_TEST_DIAGNOSTIC_FUNC*/
