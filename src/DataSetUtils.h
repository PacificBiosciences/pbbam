#ifndef PBBAM_DATASETUTILS_H
#define PBBAM_DATASETUTILS_H

#include <pbbam/Config.h>

#include <pbbam/DataSetTypes.h>

namespace PacBio {
namespace BAM {
namespace internal {

static const std::string XML_VERSION = std::string{"3.0.1"};

template <typename T>
const T& NullObject()
{
    static const T empty;
    return empty;
}

template <>
inline const PacBio::BAM::DataSetMetadata& NullObject()
{
    static const PacBio::BAM::DataSetMetadata empty("", "");
    return empty;
}

std::string GenerateUuid();

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#ifndef FETCH_CHILD_CONST_REF
#define FETCH_CHILD_CONST_REF(Class, Type, Method)            \
                                                              \
    const PacBio::BAM::Type& Class::Method() const            \
    {                                                         \
        try {                                                 \
            return Child<PacBio::BAM::Type>(#Type);           \
        } catch (std::exception&) {                           \
            return internal::NullObject<PacBio::BAM::Type>(); \
        }                                                     \
    }
#endif

#ifndef FETCH_CHILD_REF
#define FETCH_CHILD_REF(Class, Type, Method)                                       \
                                                                                   \
    PacBio::BAM::Type& Class::Method()                                             \
    {                                                                              \
        if (!HasChild(#Type)) AddChild(internal::NullObject<PacBio::BAM::Type>()); \
        return Child<PacBio::BAM::Type>(#Type);                                    \
    }
#endif

#ifndef DEFINE_ACCESSORS
#define DEFINE_ACCESSORS(Class, Type, Method)  \
    FETCH_CHILD_CONST_REF(Class, Type, Method) \
    FETCH_CHILD_REF(Class, Type, Method)
#endif

#endif  // PBBAM_DATASETUTILS_H
