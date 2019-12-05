


#include "tests/diagnostic/include.h"

using namespace PHARE::core;
using namespace PHARE::amr;
using namespace PHARE::solver;

constexpr size_t nbrGhosts = 5; // TODO UNDO REDO

struct PatchInfo
{
    std::vector<std::vector<float>> values;
    std::vector<std::string> origins;
};

template<size_t dim>
struct CoordPatchInfo
{
    using Point = PHARE::core::Point<double, dim>;
    CoordPatchInfo(std::string csv, std::vector<float> const& _vec)
        : point{Point::fromString(csv)}
        , vec{_vec}
    {
    }

    Point point;
    std::vector<float> const& vec;
};

template<size_t dim>
bool operator<(std::shared_ptr<CoordPatchInfo<dim>> const& a,
               std::shared_ptr<CoordPatchInfo<dim>> const& b)
{
    return a->point < b->point;
}

TYPED_TEST(SimulatorTest, verifyCoarsestPeriodicityOfFields)
{
    assert(TypeParam::dimension == 1); // only 1d supported for now

    TypeParam sim;

    auto& db    = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto& simdb = db["simulation"];
    EXPECT_EQ(simdb["boundary_types"].template to<std::string>(), "periodic");

    auto& hybridModel = *sim.getHybridModel();
    auto& rm          = *hybridModel.resourcesManager;
    {
        Hi5Diagnostic<TypeParam> hi5{sim, "electromag", NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.electromag("/EM_B")).addDiagDict(hi5.electromag("/EM_E")).dump();
    }
    Hi5Diagnostic<TypeParam> hi5{sim, "electromag", HighFive::File::ReadOnly};
    auto& hifile = hi5.writer.file();

    std::unordered_map<std::string, PatchInfo> patches;

    auto loadEB = [&](auto& path, std::string eb, std::string xyz) {
        std::string ebxyz   = eb + "/" + xyz;
        std::string dataset = path + ebxyz;
        if (!patches.count(ebxyz))
            patches.emplace(ebxyz, PatchInfo{});
        patches[ebxyz].values.emplace_back();
        hifile.getDataSet(dataset).read(patches[ebxyz].values.back());
        patches[ebxyz].origins.emplace_back();
        hifile.getGroup(path).getAttribute("origin").read(patches[ebxyz].origins.back());
    };

    for (auto const& time : hifile.getGroup("/").listObjectNames())
    {
        for (auto const& leaf : hifile.getGroup("/" + time + "/pl0").listObjectNames())
        {
            std::string path = "/" + time + "/pl0/" + leaf + "/";
            for (auto const& xyz : {"x", "y", "z"})
            {
                loadEB(path, "EM_E", xyz);
                loadEB(path, "EM_B", xyz);
            }
        }
    }

    std::unordered_map<std::string, std::vector<float>> values;
    for (auto const& [path, pInfo] : patches)
    {
        std::vector<std::shared_ptr<CoordPatchInfo<TypeParam::dimension>>> points;
        for (size_t i = 0; i < pInfo.origins.size(); i++)
            points.emplace_back(std::make_shared<CoordPatchInfo<TypeParam::dimension>>(
                pInfo.origins[i], pInfo.values[i]));
        std::sort(points.begin(), points.end());
        values.emplace(path, std::vector<float>{});
        for (auto const& p : points)
            for (auto const& v : p->vec)
                values[path].emplace_back(v);
    }

    auto checkDual = [&](std::vector<std::string> keys) {
        for (const auto& key : keys)
        {
            auto& v = values[key];
            for (size_t i = 0; i < nbrGhosts * 2; i++)
            {
                size_t mirror = v.size() - (nbrGhosts * 2) + i;
                EXPECT_EQ(v[i], v[mirror]);
            }
        }
    };

    auto checkPrimal = [&](std::vector<std::string> keys) {
        for (const auto& key : keys)
        {
            auto& v = values[key];
            for (size_t i = 0; i < nbrGhosts * 2 + 1; i++)
            {
                size_t mirror = v.size() - (nbrGhosts * 2) + i - 1;
                EXPECT_EQ(v[i], v[mirror]);
            }
        }
    };

    checkDual({"EM_E/x", "EM_B/y", "EM_B/z"});
    checkPrimal({"EM_B/x", "EM_E/y", "EM_E/z"});
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
