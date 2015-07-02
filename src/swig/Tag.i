/* Tag.i */

%module PacBioBam

%{
#include <pbbam/Tag.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::Tag::Tag(Tag&&);
%ignore PacBio::BAM::Tag::operator=;

HANDLE_STD_EXCEPTION(ToInt8);
HANDLE_STD_EXCEPTION(ToUInt8);
HANDLE_STD_EXCEPTION(ToInt16);
HANDLE_STD_EXCEPTION(ToUInt16);
HANDLE_STD_EXCEPTION(ToInt32);
HANDLE_STD_EXCEPTION(ToUInt32);
HANDLE_STD_EXCEPTION(ToFloat);
HANDLE_STD_EXCEPTION(ToString);
HANDLE_STD_EXCEPTION(ToInt8Array);
HANDLE_STD_EXCEPTION(ToUInt8Array);
HANDLE_STD_EXCEPTION(ToInt16Array);
HANDLE_STD_EXCEPTION(ToUInt16Array);
HANDLE_STD_EXCEPTION(ToInt32Array);
HANDLE_STD_EXCEPTION(ToUInt32Array);
HANDLE_STD_EXCEPTION(ToFloatArray);

#ifdef SWIGR

%ignore PacBio::BAM::Tag::Tag(int8_t value);
%ignore PacBio::BAM::Tag::Tag(uint8_t value);
%ignore PacBio::BAM::Tag::Tag(int16_t value);
%ignore PacBio::BAM::Tag::Tag(uint16_t value);
%ignore PacBio::BAM::Tag::Tag(int32_t value);
%ignore PacBio::BAM::Tag::Tag(uint32_t value);
%ignore PacBio::BAM::Tag::Tag(float value);
%ignore PacBio::BAM::Tag::Tag(const std::string& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<int8_t>& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<uint8_t>& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<int16_t>& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<uint16_t>& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<int32_t>& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<uint32_t>& value);
%ignore PacBio::BAM::Tag::Tag(const std::vector<float>& value);

%extend PacBio::BAM::Tag {
	
	PacBio::BAM::Tag FromInt8(int x)   { return PacBio::BAM::Tag(static_cast<int8_t>(x));   }
	PacBio::BAM::Tag FromUInt8(int x)  { return PacBio::BAM::Tag(static_cast<uint8_t>(x));  }
	PacBio::BAM::Tag FromInt16(int x)  { return PacBio::BAM::Tag(static_cast<int16_t>(x));  }
	PacBio::BAM::Tag FromUInt16(int x) { return PacBio::BAM::Tag(static_cast<uint16_t>(x)); }
	PacBio::BAM::Tag FromInt32(int x)  { return PacBio::BAM::Tag(static_cast<int32_t>(x));  }
	PacBio::BAM::Tag FromUInt32(int x) { return PacBio::BAM::Tag(static_cast<uint32_t>(x)); }
	PacBio::BAM::Tag FromFloat(int x)  { return PacBio::BAM::Tag(static_cast<float>(x));    }
	
	PacBio::BAM::Tag FromInt8Array(const std::vector<int>& v)
	{
		std::vector<int8_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<int8_t>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
	
	PacBio::BAM::Tag FromUInt8Array(const std::vector<int>& v)
	{
		std::vector<uint8_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<uint8_t>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
	
	PacBio::BAM::Tag FromInt16Array(const std::vector<int>& v)
	{
		std::vector<int16_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<int16_t>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
	
	PacBio::BAM::Tag FromUInt16Array(const std::vector<int>& v)
	{
		std::vector<int16_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<uint16_t>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
	
	PacBio::BAM::Tag FromInt32Array(const std::vector<int>& v)
	{
		std::vector<int16_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<int32_t>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
	
	PacBio::BAM::Tag FromUInt32Array(const std::vector<int>& v)
	{
		std::vector<int16_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<uint32_t>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
	
	PacBio::BAM::Tag FromFloatArray(const std::vector<int>& v)
	{
		std::vector<int16_t> result;
		const size_t numElements = v.size();
		result.reserve(numElements);
		for (size_t i = 0; i < numElements; ++i) 
			result.push_back(static_cast<float>(v.at(i)));
		return PacBio::BAM::Tag(result); 
	}
}
#endif // SWIGR

%include <pbbam/Tag.h>