
template<std::size_t dim>
using NdArray_t = typename PHARE::core::PHARE_Types<dim>::Array_t;
using NdArrays  = ::testing::Types<NdArray_t<1>, NdArray_t<2>, NdArray_t<3>>;

template<typename NdArray>
struct NdArrayTest : public ::testing::Test
{
    static const auto dimension = NdArray::dimension;
};
TYPED_TEST_SUITE(NdArrayTest, NdArrays);

TYPED_TEST(NdArrayTest, is_raja_exec_settable)
{
    static constexpr auto dimension = TypeParam::dimension;
    static constexpr auto shape     = ConstArray<std::uint32_t, dimension>(10);

    TypeParam src_arr{shape}, dst_arr{shape};
    for (auto& e : src_arr)
        e = 12.;
    auto dst = dst_arr.data();
    auto src = src_arr.data();
    PHARE::raja::exec([=] RAJA_DEVICE(int i) { dst[i] = src[i]; }, dst_arr.size());
    for (auto const& e : dst_arr)
        EXPECT_DOUBLE_EQ(12., e);
}

TYPED_TEST(NdArrayTest, is_raja_copyable)
{
    static constexpr auto dimension = TypeParam::dimension;
    static constexpr auto shape     = ConstArray<std::uint32_t, dimension>(10);

    TypeParam src_arr{shape}, dst_arr{shape};
    for (auto& e : src_arr)
        e = 12.;
    PHARE::raja::copy(dst_arr.data(), src_arr.data(), dst_arr.size());
    for (auto const& e : dst_arr)
        EXPECT_DOUBLE_EQ(12., e);
}

TYPED_TEST(NdArrayTest, is_raja_settable)
{
    static constexpr auto dimension = TypeParam::dimension;
    static constexpr auto shape     = ConstArray<std::uint32_t, dimension>(10);

    TypeParam dst_arr{shape};
    PHARE::raja::set(dst_arr.data(), 12., dst_arr.size());
    for (auto const& e : dst_arr)
        EXPECT_DOUBLE_EQ(12., e);
}
