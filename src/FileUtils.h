// Author: Derek Barnett

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <chrono>
#include <string>

namespace PacBio {
namespace BAM {
namespace internal {

struct FileUtils
{
public:
    /// \returns application's current working directory
    static std::string CurrentWorkingDirectory();

    /// Parses a filepath for the the directory name for a file.
    ///
    /// Essentially this method strips the filename from the string provided (/path/to/file => /path/to).
    /// If only a filename is provided, then "." is returned to indicate the current directory.
    ///
    /// \param[in] file name of file (can be just a filename or path/to/filename)
    /// \returns file's directory name
    ///
    static std::string DirectoryName(const std::string& file);

    /// Check for existence of a file.
    ///
    /// \param[in] fn full path to file
    /// \returns true if file exists & can be opened
    ///
    static bool Exists(const char* fn);

    /// Check for existence of a file.
    ///
    /// \param[in] fn full path to file
    /// \returns true if file exists & can be opened
    ///
    static bool Exists(const std::string& fn);

    /// Check "last modified" timestamp for a file.
    ///
    /// \param[in] fn full path to file
    /// \returns time of last modification
    /// \throws runtime_error if file info can't be accessed
    ///
    static std::chrono::system_clock::time_point LastModified(const char* fn);

    /// Check "last modified" timestamp for a file.
    ///
    /// \param[in] fn full path to file
    /// \returns time of last modification
    /// \throws runtime_error if file info can't be accessed
    ///
    static std::chrono::system_clock::time_point LastModified(const std::string& fn);

    /// Resolves input file path using optional starting directory.
    ///
    /// \verbatim
    ///   /absolute/path/to/file.txt   => /absolute/path/to/file.txt
    ///   ../relative/path/to/file.txt => <from>/../relative/path/to/file.txt
    ///   file.txt                     => <from>/file.txt
    /// \endverbatim
    ///
    /// \note This method will strip any URI scheme as well ("file://") so that the result is immediately ready from I/O operations.
    ///
    /// \param[in] filePath file path to be resolved
    /// \param[in] from     optional starting directory (useful if not same as application's working directory)
    /// \returns resolved file path
    ///
    static std::string ResolvedFilePath(const std::string& filePath, const std::string& from = ".");

    /// \returns native path separator
    constexpr static char Separator();

    /// Check size of file.
    ///
    /// \param[in] fn full path to file
    /// \returns file size in bytes
    /// \throws runtime_error if file info can't be accessed
    ///
    static off_t Size(const char* fn);

    /// Check size of file.
    ///
    /// \param[in] fn full path to file
    /// \returns file size in bytes
    /// \throws runtime_error if file info can't be accessed
    ///
    static off_t Size(const std::string& fn);
};

inline bool FileUtils::Exists(const std::string& fn) { return FileUtils::Exists(fn.c_str()); }

inline std::chrono::system_clock::time_point FileUtils::LastModified(const std::string& fn)
{
    return FileUtils::LastModified(fn.c_str());
}

inline off_t FileUtils::Size(const std::string& fn) { return FileUtils::Size(fn.c_str()); }

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // FILEUTILS_H
