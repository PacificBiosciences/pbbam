#ifndef SIMPLEREADIMPL_H
#define SIMPLEREADIMPL_H

#include <cstddef>

namespace PacBio {
namespace BAM {

class SimpleRead;
class MappedSimpleRead;

namespace internal {

struct ClipResult;

void ClipSimpleRead(SimpleRead& read, const internal::ClipResult& result, size_t start, size_t end);

// NOTE: 'result' is moved into here, so we can take the CIGAR
void ClipMappedRead(MappedSimpleRead& read, internal::ClipResult result, size_t start, size_t end);

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // SIMPLEREADIMPL_H