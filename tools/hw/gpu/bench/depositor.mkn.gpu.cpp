#include <cassert>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define RNM (double)rand() / ((double)RAND_MAX)

#include "kul/gpu.hpp"
#include "kul/gpu/tuple.hpp"

static constexpr size_t X = 1024, Y = 1024, Z = 1;         // 40;
static constexpr size_t TPB_X = 16, TPB_Y = 16, TPB_Z = 1; // 4;
static constexpr uint32_t NUM        = X * Y * Z;
static constexpr uint32_t ITERATIONS = 1000;


template<typename Float_, std::size_t ppc>
__global__ void gpu_fn(Float_* x, Float_* y, Float_* dens)
{
    auto p_block_idx = kul::gpu::idx() * ppc;

    uint32_t ixy1, ixy2, ixy3, ixy4;
    uint32_t ix, iy;
    uint32_t ixd, iyd;

    uint32_t repeat = ITERATIONS;

    Float_ dix, diy;
    Float_ odx, ody;
    Float_ w1, w2, w3, w4;

    do
    {
        for (uint32_t ip = p_block_idx; ip < p_block_idx + ppc; ip++)
        {
            ixd = x[ip] * odx + 1;
            iyd = y[ip] * ody + 1;

            ix = (uint32_t)ixd;
            iy = (uint32_t)iyd;

            dix = ixd - ix;
            diy = iyd - iy;

            w1 = (1.0 - dix) * (1.0 - diy);
            w2 = (1.0 - dix) * (diy);
            w3 = (dix) * (diy);
            w4 = (dix) * (1.0 - diy);

            ixy1 = iy + (ix)*Y;
            ixy2 = iy + 1 + (ix)*Y;
            ixy3 = iy + 1 + (ix + 1) * Y;
            ixy4 = iy + (ix + 1) * Y;

            dens[ixy1] += w1;
            dens[ixy2] += w2;
            dens[ixy3] += w3;
            dens[ixy4] += w4;
        }
        --repeat;
    } while (repeat != 0);
}

template<typename Float_, std::size_t ppc>
void do_thing()
{
    constexpr Float_ xm = 40, ym = 40;
    constexpr uint32_t npart = X * Y * ppc;

    std::vector<Float_> dens(NUM, 1), x(npart), y(npart);

    for (uint32_t ip = 0; ip < npart; ip++)
    {
        x[ip] = RNM * xm;
        y[ip] = RNM * ym;
    }

    clock_t start, end;

    kul::gpu::DeviceMem<Float_> devDens(dens), devX(x), devY(y);
    start = clock();
    kul::gpu::Launcher{X, Y, Z, TPB_X, TPB_Y, TPB_Z}(gpu_fn<Float_, ppc>, devX, devY, devDens);
    end = clock();

    auto sec = (float)(end - start) / CLOCKS_PER_SEC;

    auto h_dens = devDens();

    assert(h_dens.size() == NUM);
    assert(h_dens != dens);

    printf("time = %f seconds with %lu ppc\n", sec, ppc);
}


int main(int argc, char** argv)
{
    std::cout << "float32 launching" << std::endl;
    do_thing<float, 50>();
    do_thing<float, 100>();
    do_thing<float, 128>();
    do_thing<float, 256>();
    do_thing<float, 500>();
    do_thing<float, 1000>();
    std::cout << "float32 finished" << std::endl;

    std::cout << "float64 launching" << std::endl;
    do_thing<double, 50>();
    do_thing<double, 100>();
    do_thing<double, 128>();
    do_thing<double, 256>();
    do_thing<double, 500>();
    std::cout << "float64 finished" << std::endl;

    return 0;
}
