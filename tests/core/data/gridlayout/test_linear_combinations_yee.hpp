#ifndef PHARE_TEST_LINEAR_COMBINAISON_HPP
#define PHARE_TEST_LINEAR_COMBINAISON_HPP


#include <array>
#include <vector>
#include <cstddef>
#include <fstream>



struct Combination
{
    int dimension;
    int interpOrder;
    std::vector<int> ix, iy, iz;
    double coef;
};



std::array<Combination, 9> inline readFile(std::string const& filename)
{
    std::ifstream file{filename};
    std::array<Combination, 9> combinations;
    int nbrPts;

    for (std::size_t c = 0; c < 9; ++c)
    {
        file >> combinations[c].dimension >> combinations[c].interpOrder >> nbrPts;

        for (int ip = 0; ip < nbrPts; ++ip)
        {
            int i, j, k;
            if (combinations[c].dimension == 1)
            {
                file >> i;
                combinations[c].ix.push_back(i);
            }
            else if (combinations[c].dimension == 2)
            {
                file >> i >> j;
                combinations[c].ix.push_back(i);
                combinations[c].iy.push_back(j);
            }
            else
            {
                file >> i >> j >> k;
                combinations[c].ix.push_back(i);
                combinations[c].iy.push_back(j);
                combinations[c].iz.push_back(k);
            }
        }
        file >> combinations[c].coef;
    }


    return combinations;
}




#endif
