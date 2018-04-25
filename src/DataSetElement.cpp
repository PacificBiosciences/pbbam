// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/internal/DataSetElement.h"

#include "DataSetUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

const std::string& DataSetElement::SharedNullString()
{
    return internal::NullObject<std::string>();
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
