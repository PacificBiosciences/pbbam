// File Description
/// \file PbiRawData.inl
/// \brief Inline implementations for the classes used for working with raw PBI
///        data.
//
// Author: Derek Barnett

#include "pbbam/PbiRawData.h"

namespace PacBio {
namespace BAM {

inline const PbiRawBarcodeData& PbiRawData::BarcodeData() const
{ return barcodeData_; }

inline PbiRawBarcodeData& PbiRawData::BarcodeData()
{ return barcodeData_; }

inline const PbiRawBasicData& PbiRawData::BasicData() const
{ return basicData_; }

inline PbiRawBasicData& PbiRawData::BasicData()
{ return basicData_; }

inline std::string PbiRawData::Filename() const
{ return filename_; }

inline PbiFile::Sections PbiRawData::FileSections() const
{ return sections_; }

inline PbiRawData& PbiRawData::FileSections(PbiFile::Sections sections)
{ sections_ = sections; return *this; }

inline bool PbiRawData::HasBarcodeData() const
{ return HasSection(PbiFile::BARCODE); }

inline bool PbiRawData::HasMappedData() const
{ return HasSection(PbiFile::MAPPED); }

inline bool PbiRawData::HasReferenceData() const
{ return HasSection(PbiFile::REFERENCE); }

inline bool PbiRawData::HasSection(const PbiFile::Section section) const
{ return (sections_ & section) != 0; }

inline uint32_t PbiRawData::NumReads() const
{ return numReads_; }

inline PbiRawData& PbiRawData::NumReads(uint32_t num)
{ numReads_ = num; return *this; }

inline const PbiRawMappedData& PbiRawData::MappedData() const
{ return mappedData_; }

inline PbiRawMappedData& PbiRawData::MappedData()
{ return mappedData_; }

inline const PbiRawReferenceData& PbiRawData::ReferenceData() const
{ return referenceData_; }

inline PbiRawReferenceData& PbiRawData::ReferenceData()
{ return referenceData_; }

inline PbiFile::VersionEnum PbiRawData::Version() const
{ return version_; }

inline PbiRawData& PbiRawData::Version(PbiFile::VersionEnum version)
{ version_ = version; return *this; }

inline bool PbiReferenceEntry::operator==(const PbiReferenceEntry& other) const
{
    return tId_      == other.tId_ &&
           beginRow_ == other.beginRow_ &&
           endRow_   == other.endRow_;
}

} // namespace BAM
} // namespace PacBio
