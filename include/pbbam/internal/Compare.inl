#ifndef PBBAM_COMPARE_INL
#define PBBAM_COMPARE_INL

#include <pbbam/Config.h>

#include <pbbam/Compare.h>

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T, T>
struct MemberFnProxy;

template <typename T, typename R, typename... Args, R (T::*fn)(Args...) const>
struct MemberFnProxy<R (T::*)(Args...) const, fn>
{
    static R call(const T& obj, Args&&... args) { return (obj.*fn)(std::forward<Args>(args)...); }
};

}  // namespace internal

template <typename ValueType,
          typename Compare::MemberFunctionBaseHelper<ValueType>::MemberFnType fn,
          typename CompareType>
bool Compare::MemberFunctionBase<ValueType, fn, CompareType>::operator()(const BamRecord& lhs,
                                                                         const BamRecord& rhs) const
{
    using MemberFnTypeImpl = typename Compare::MemberFunctionBaseHelper<ValueType>::MemberFnType;
    using Proxy = internal::MemberFnProxy<MemberFnTypeImpl, fn>;

    CompareType cmp;
    return cmp(Proxy::call(lhs), Proxy::call(rhs));
}

inline bool Compare::None::operator()(const BamRecord&, const BamRecord&) const noexcept
{
    return false;
}

}  // namespace BAM
}  // namespace PacBio

#endif // PBBAM_COMPARE_INL
