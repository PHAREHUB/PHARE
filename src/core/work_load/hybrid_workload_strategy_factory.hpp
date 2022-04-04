

class HybridWorkLoadStrategyFactory
{
    public :
        static std::unique_ptr<HybridWorkLoadEstimatorStrategy> create(std::string stratName)
            {
                if (stratName == "NPPC")
                    return std::make_unique<ConcreteHybridWorkLoadEstimatorStrategyNPPC>();
                else
                    return {};
            };
};

