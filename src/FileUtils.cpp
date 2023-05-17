#include "PbbamInternalConfig.h"

#include "FileUtils.h"

#include <pbbam/StringUtilities.h>

#include <boost/algorithm/string.hpp>

#include <memory>
#include <sstream>
#include <stdexcept>

#include <cassert>
#include <cstddef>

#include <sys/stat.h>
#include <unistd.h>

namespace PacBio {
namespace BAM {

namespace {

// pops "file://" scheme off the front of a URI/filepath, if found
static std::string removeFileUriScheme(const std::string& uri)
{
    assert(!uri.empty());

    auto schemeLess = uri;
    const auto fileScheme = std::string{"file://"};
    const auto schemeFound = schemeLess.find(fileScheme);
    if (schemeFound != std::string::npos) {
        if (schemeFound != 0) {
            std::ostringstream s;
            s << "[pbbam] file utilities ERROR: malformed URI, scheme is not at beginning:\n"
              << "  uri: " << uri;
            throw std::runtime_error{s.str()};
        }
        schemeLess = schemeLess.substr(fileScheme.size());
    }
    return schemeLess;
}

#ifdef PBBAM_WIN_FILEPATHS

static std::string removeDiskName(const std::string& filePath)
{
    if (filePath.size() >= 2) {
        const char firstChar = filePath.at(0);
        if ((isalpha(firstChar) != 0) && (filePath.at(1) == ':')) return filePath.substr(2);
    }
    return filePath;
}

constexpr char NATIVE_PATH_SEPARATOR = '\\';

static bool native_pathIsAbsolute(const std::string& filePath)
{
    assert(!filePath.empty());

    // if starts with single slash or double slash
    if (boost::algorithm::starts_with(filePath, "\\")) return true;

    // if starts with single or double-dots -> not absolute
    if (boost::algorithm::starts_with(filePath, ".")) return false;

    // if starts with disk drive name and colon ("C:\foo\bar.txt")
    // strip the drive name and check to see if the remaining path is absolute
    if (filePath.size() >= 2) {
        const char firstChar = filePath.at(0);
        if ((isalpha(firstChar) != 0) && (filePath.at(1) == ':'))
            return native_pathIsAbsolute(removeDiskName(filePath));
    }

    // otherwise, likely relative
    return false;
}

static std::string native_resolvedFilePath(const std::string& filePath, const std::string& from)
{
    // strip file:// scheme if present
    auto schemeLess = removeFileUriScheme(filePath);

    // if empty or already absolute path, just return it
    // upfront empty check simplifies further parsing logic
    if (schemeLess.empty() || native_pathIsAbsolute(schemeLess)) return schemeLess;

    // else make relative from the provided 'from' directory
    //
    // first pop disk name, then any leading single-dot '.'
    //
    // since we're prepending the 'from' directory, we can remove
    // any leading './' form our file path. this may just mean that
    // we pop it off to add it right back (when from == '.'), but this
    // keeps it consistent with other 'from' parent directories
    //
    schemeLess = removeDiskName(schemeLess);

    const bool thisDirAtStart = (schemeLess.find(".") == 0);
    if (thisDirAtStart) {
        if (schemeLess.find(NATIVE_PATH_SEPARATOR) == 1) schemeLess = schemeLess.substr(2);
    }
    return from + NATIVE_PATH_SEPARATOR + schemeLess;
}

#else  // else for non-Windows systems

constexpr char NATIVE_PATH_SEPARATOR = '/';

static bool native_pathIsAbsolute(const std::string& filePath) { return filePath.at(0) == '/'; }

static std::string native_resolvedFilePath(const std::string& filePath, const std::string& from)
{
    // strip file:// scheme if present
    auto schemeLess = removeFileUriScheme(filePath);

    // if empty or already absolute path, just return it
    // upfront empty check simplifies further parsing logic
    if (schemeLess.empty() || native_pathIsAbsolute(schemeLess)) {
        return schemeLess;
    }

    // else make relative from the provided 'from' directory
    //
    // since we're prepending the 'from' directory, we can remove
    // any leading './' form our file path. this may just mean that
    // we pop it off to add it right back (when from == '.'), but this
    // keeps it consistent with other 'from' parent directories
    //
    const bool thisDirAtStart = (schemeLess.find(".") == 0);
    if (thisDirAtStart) {
        if (schemeLess.find(NATIVE_PATH_SEPARATOR) == 1) {
            schemeLess = schemeLess.substr(2);
        }
    }
    return from + NATIVE_PATH_SEPARATOR + schemeLess;
}

#endif  // PBBAM_WIN_FILEPATHS

}  // namespace

// see http://stackoverflow.com/questions/2869594/how-return-a-stdstring-from-cs-getcwd-function
std::string FileUtils::CurrentWorkingDirectory()
{
    const std::size_t chunkSize = 1024;
    const std::size_t maxNumChunks = 20;

    // stack-based buffer for 'normal' case
    char buffer[chunkSize];
    if (getcwd(buffer, sizeof(buffer)) != nullptr) {
        return std::string(buffer);
    }

    // if error is not ERANGE, then it's not a problem of too-long name... something else happened
    if (errno != ERANGE) {
        throw std::runtime_error{
            "[pbbam] file utilities ERROR: could not determine current working directory path"};
    }

    // long path - use heap, trying progressively longer buffers
    for (std::size_t chunks = 2; chunks < maxNumChunks; ++chunks) {
        auto cwd = std::make_unique<char[]>(chunkSize * chunks);
        if (getcwd(cwd.get(), chunkSize * chunks) != nullptr) {
            return std::string(cwd.get());
        }

        // if error is not ERANGE, then it's not a problem of too-long name... something else happened
        if (errno != ERANGE) {
            throw std::runtime_error{
                "[pbbam] file utilities ERROR: could not determine current working directory path"};
        }
    }

    // crazy long path name
    throw std::runtime_error{
        "[pbbam] file utilities ERROR: could not determine current working directory - extremely "
        "long path"};
}

std::string FileUtils::DirectoryName(const std::string& file)
{
    const auto found = file.rfind(Separator(), file.length());
    if (found != std::string::npos) {
        return file.substr(0, found);
    }
    return std::string(".");
}

bool FileUtils::Exists(const char* fn)
{
    struct stat buf;
    return (stat(fn, &buf) != -1);
}

std::chrono::system_clock::time_point FileUtils::LastModified(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0) {
        std::ostringstream msg;
        msg << "[pbbam] file utilities ERROR: could not determine 'last modified' timestamp:\n"
            << "  file: " << fn;
        throw std::runtime_error{msg.str()};
    }
    return std::chrono::system_clock::from_time_t(s.st_mtime);
}

std::string FileUtils::ResolvedFilePath(const std::string& filePath, const std::string& from)
{
    return native_resolvedFilePath(filePath, from);
}

char FileUtils::Separator() { return NATIVE_PATH_SEPARATOR; }

off_t FileUtils::Size(const char* fn)
{
    struct stat s;
    if (stat(fn, &s) != 0) {
        std::ostringstream msg;
        msg << "[pbbam] file utilities ERROR: could not determine file size:\n"
            << "  file: " << fn;
        throw std::runtime_error{msg.str()};
    }
    return s.st_size;
}

}  // namespace BAM
}  // namespace PacBio
