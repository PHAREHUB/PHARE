
#ifndef PHARE_TEST_DIAGNOSTIC_FUNC
#define PHARE_TEST_DIAGNOSTIC_FUNC

namespace PHARE_test
{
double density(double x)
{
    return x * x + 2.;
}
std::array<double, 3> bulkVelocity([[maybe_unused]] double x)
{
    return std::array<double, 3>{{1.0, 0.0, 0.0}};
}
std::array<double, 3> thermalVelocity([[maybe_unused]] double x)
{
    return std::array<double, 3>{{0.5, 0.0, 0.0}};
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