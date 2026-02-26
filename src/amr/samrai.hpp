#ifndef PHARE_AMR_SAMRAI_HPP
#define PHARE_AMR_SAMRAI_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"

#include <SAMRAI/tbox/RestartManager.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include <SAMRAI/hier/PatchDataRestartManager.h>

#include <iostream>

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(std::string const& message, std::string const& filename, int const line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle //
{
public:
    SamraiLifeCycle(int argc = 0, char** argv = nullptr);

    ~SamraiLifeCycle();

    static void reset();

    static SAMRAI::hier::VariableDatabase* getDatabase();

    static SAMRAI::hier::PatchDataRestartManager* getPatchDataRestartManager();

    static SAMRAI::tbox::RestartManager* getRestartManager();
};


} // namespace PHARE

namespace PHARE::amr
{


template<typename T, typename A>
auto& getVectorFromRestart(auto& db, auto const& path, std::vector<T, A>& vec)
{
    if (vec.size() == 0)
        throw std::runtime_error("SAMRAI Restarts: vectors must be presized as expected");

    if constexpr (std::is_same_v<T, double>)
        db.getDoubleArray(path, vec.data(), vec.size());

    else if constexpr (std::is_same_v<T, int>)
        db.getIntegerArray(path, vec.data(), vec.size());

    else
        static_assert(core::dependent_false_v<T>,
                      "SAMRAI getFromRestart Vector not supported, add it!");

    return vec;
};

template<typename T, typename A>
void putVectorToRestart(auto& db, auto const& path, std::vector<T, A> const& vec)
{
    if constexpr (std::is_same_v<T, double>)
        db.putDoubleArray(path, vec.data(), vec.size());

    else if constexpr (std::is_same_v<T, int>)
        db.putIntegerArray(path, vec.data(), vec.size());

    else
        static_assert(core::dependent_false_v<T>,
                      "SAMRAI putToRestart Vector not supported, add it!");
};


} // namespace PHARE::amr

#endif /*PHARE_AMR_SAMRAI_HPP*/
