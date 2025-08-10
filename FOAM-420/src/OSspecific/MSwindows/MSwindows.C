/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2016 blueCAPE Ltda.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2011 Symscape
    (c) 2011-2017 FSD blueCAPE Lda

Description
    MS Windows versions of the functions declared in OSspecific.H

Details
    This file has been created by blueCAPE's unofficial mingw patches for
    OpenFOAM.
    For more information about these patches, visit:
        http://bluecfd.com/Core

    Details on how this file was created:
      - This file was originally based on Symscape's own work on patching
        OpenFOAM for working on Windows, circa 2009.
      - Further changes were made by blueCAPE over the years.
      - Includes code from POSIX.C where applicable.


\*---------------------------------------------------------------------------*/

#include "include/OSspecific.H"
#include "MSwindows.H"
#include "global/foamVersion.H"
#include "primitives/strings/fileName/fileName.H"
#include "fileStat.H"
#include "timer.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOstreams/Pstreams/Pstream.H"

// Undefine DebugInformation, because we don't need it and it collides with a macro
// in windows.h
#undef DebugInformation

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <map>

// Windows system header files
#include <io.h> // _close
#include <windows.h>
#include <signal.h>

// For C++11 std::thread and std::mutex
#include <thread>
#include <mutex>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(MSwindows, 0);
}

namespace Foam
{

// Don't abort under windows, causes abort dialog to
// popup. Instead just exit with exitCode.
static
void sigAbortHandler(int exitCode)
{
  ::exit(exitCode);
}


static
bool installAbortHandler()
{
  // If it didn't succeed there's not much we can do,
  // so don't check result.
  ::signal(SIGABRT, &sigAbortHandler);
  return true;
}


static bool const abortHandlerInstalled = installAbortHandler();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Get last windows api error from GetLastError
string MSwindows::getLastError()
{
    // Based on an example at:
    // http://msdn2.microsoft.com/en-us/library/ms680582(VS.85).aspx

    LPVOID lpMsgBuf;
    LPVOID lpDisplayBuf;
    DWORD dw = GetLastError();

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        nullptr,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        LPTSTR(&lpMsgBuf),
        0, nullptr );

    lpDisplayBuf = LocalAlloc(LMEM_ZEROINIT,
        (lstrlen(static_cast<LPCTSTR>(lpMsgBuf))+40)*sizeof(TCHAR));
    sprintf(static_cast<LPTSTR>(lpDisplayBuf),
            "Error %d: %s", int(dw), static_cast<LPCTSTR>(lpMsgBuf));

    const string errorMessage = static_cast<LPTSTR>(lpDisplayBuf);

    LocalFree(lpMsgBuf);
    LocalFree(lpDisplayBuf);

    return errorMessage;
}

} // namespace Foam

namespace Foam
{
//-Declared here to avoid polluting MSwindows.H with windows.h
namespace MSwindows
{
    //- Get windows user name
    Foam::string getUserName();

    //- Remove quotes, if any, from string
    void removeQuotes(Foam::string & arg);

    //- Convert windows directory slash (back-slash) to unix (forward-slash).
    //- Windows is fine with unix like directory slashes.
    //- Foam's file io (see src/OpenFOAM/db/IOstreams/Sstreams/OSwrite.C)
    //- uses back-slash as escape character and continuation,
    //- so not an option to have windows file paths with back-slashes
    void toUnixSlash(Foam::string & arg);

    //- Auto create and then delete array when this goes out of scope
    template<class T>
    class AutoArray
    {
      T* const array_;

    public:
      AutoArray(const unsigned long arrayLength);
      ~AutoArray();

      //- Access array
      T* get();
    }; // class AutoArray


    //- Directory contents iterator
    class DirectoryIterator
    {
      WIN32_FIND_DATA findData_;
      HANDLE findHandle_;
      fileName nextName_;
      bool hasMore_;

    public:
      DirectoryIterator(const fileName & directory);
      ~DirectoryIterator();

      //- Initialization succeeded
      bool isValid() const;

      //- Has more?
      bool hasNext() const;

      //- Next item
      const fileName & next();
    }; // class DirectoryIterator
} // namespace MSwindows
} // namespace Foam

inline
void Foam::MSwindows::removeQuotes(Foam::string & arg)
{
    std::size_t pos;

    while (Foam::string::npos != (pos = arg.find('"')))
    {
        arg.erase(pos, 1);
    }
}


inline
void Foam::MSwindows::toUnixSlash(Foam::string & arg)
{
    arg.replaceAll("\\", "/");

    const string UNC("//");

    // Preserve UNC i.e., \\machine-name\...
    if (0 == arg.find(UNC))
    {
        arg.replace(UNC, "\\\\");
    }
}


Foam::string Foam::MSwindows::getUserName()
{
    const DWORD bufferSize = 256;
    TCHAR buffer[bufferSize];
    DWORD actualBufferSize = bufferSize;
    string nameAsString;

    bool success = ::GetUserName(buffer, &actualBufferSize);

    if (success)
    {
        nameAsString = buffer;
    }
    else
    {
        if (ERROR_INSUFFICIENT_BUFFER == ::GetLastError() &&
            32768 > actualBufferSize)
        {
            AutoArray<TCHAR> actualBuffer(actualBufferSize);
            ::GetUserName(actualBuffer.get(), &actualBufferSize);
            nameAsString = actualBuffer.get();
        }
    }

    return nameAsString;
}


template<class T>
inline
Foam::MSwindows::AutoArray<T>::AutoArray(const unsigned long arrayLength)
    : array_(new T[arrayLength])
{}


template<class T>
inline
Foam::MSwindows::AutoArray<T>::~AutoArray()
{
    delete [] array_;
}


template<class T>
inline
T* Foam::MSwindows::AutoArray<T>::get()
{
    return array_;
}


inline
bool Foam::MSwindows::DirectoryIterator::isValid() const
{
    const bool valid = (INVALID_HANDLE_VALUE != findHandle_);
    return valid;
}


Foam::MSwindows::DirectoryIterator::DirectoryIterator(const fileName & directory)
{
    const fileName directoryContents = directory/"*";
    findHandle_ = ::FindFirstFile(directoryContents.c_str(), &findData_);
    hasMore_    = isValid();
}


Foam::MSwindows::DirectoryIterator::~DirectoryIterator()
{
    if (isValid())
    {
        ::FindClose(findHandle_);
    }
}


inline
bool Foam::MSwindows::DirectoryIterator::hasNext() const
{
    assert(isValid());

    return hasMore_;
}


inline
const Foam::fileName & Foam::MSwindows::DirectoryIterator::next()
{
    assert(hasNext());

    nextName_ = findData_.cFileName;
    hasMore_  = ::FindNextFile(findHandle_, &findData_);

    return nextName_;
}


pid_t Foam::pid()
{
#ifdef WIN32
    const DWORD processId = ::GetCurrentProcessId();
#elif WIN64
    const pid_t processId = (pid_t) ::GetCurrentProcessId();
#endif
    return processId;
}


pid_t Foam::ppid()
{
    // No equivalent under windows.

    if (MSwindows::debug)
    {
        InfoInFunction
            << "ppid not supported under MSwindows" << endl;
    }

    return 0;
}


pid_t Foam::pgid()
{
    // No equivalent under windows.

    if (MSwindows::debug)
    {
        InfoInFunction
            << "pgid not supported under MSwindows" << endl;
    }

    return 0;
}


bool Foam::env(const word& envName)
{
    const DWORD actualBufferSize =
      ::GetEnvironmentVariable(envName.c_str(), nullptr, 0);

    const bool envExists = (0 < actualBufferSize);
    return envExists;
}


Foam::string Foam::getEnv(const word& envName)
{
    string envAsString;

    const DWORD actualBufferSize =
      ::GetEnvironmentVariable(envName.c_str(), nullptr, 0);

    if (0 < actualBufferSize)
    {
        MSwindows::AutoArray<TCHAR> actualBuffer(actualBufferSize);
        ::GetEnvironmentVariable(envName.c_str(),
                                 actualBuffer.get(),
                                 actualBufferSize);
        envAsString = actualBuffer.get();
        toUnixPath(envAsString);
    }

    return envAsString;
}


bool Foam::setEnv
(
    const word& envName,
    const std::string& value,
    const bool /*overwrite*/
)
{
    const bool success =
      ::SetEnvironmentVariable(envName.c_str(), value.c_str());
    return success;
}


Foam::string Foam::hostName(bool full)
{
    const DWORD bufferSize = MAX_COMPUTERNAME_LENGTH + 1;
    TCHAR buffer[bufferSize];
    DWORD actualBufferSize = bufferSize;

    const bool success =
      ::GetComputerName(buffer, &actualBufferSize);
    const string computerName = success ? buffer : string::null;
    return computerName;
}


Foam::string Foam::domainName()
{
    // FIXME: this should be implemented for completion.
    // Could use ::gethostname and ::gethostbyname like POSIX.C, but would
    // then need to link against ws_32. Prefer to minimize dependencies.
    // See task https://github.com/blueCFD/Core/issues/61

    return string::null;
}


Foam::string Foam::userName()
{
    string nameAsString = getEnv((word)"USERNAME");

    if (nameAsString.empty())
    {
        nameAsString = MSwindows::getUserName();
    }

    return nameAsString;
}


bool Foam::isAdministrator()
{
    // FIXME: This is going to be needed as well... for building dynamic code.
    // Going to assume to be false, so this way it's possible to build dynamic
    // code!
    // Besides, in Windows XP the usual default is to be the administrator...
    // See task https://github.com/blueCFD/Core/issues/62

    return false;
}

Foam::fileName Foam::home()
{
    string homeDir = getEnv((word)"HOME");

    if (homeDir.empty())
    {
        homeDir = getEnv((word)"USERPROFILE");
    }

    return homeDir;
}


Foam::fileName Foam::home(const std::string& userName)
{
    return home();
}


Foam::fileName Foam::cwd()
{
    string currentDirectory;

    const DWORD actualBufferSize =
      ::GetCurrentDirectory(0, nullptr);

    if (0 < actualBufferSize)
    {
        MSwindows::AutoArray<TCHAR> actualBuffer(actualBufferSize);
        ::GetCurrentDirectory(actualBufferSize,
                              actualBuffer.get());
        currentDirectory = actualBuffer.get();
        MSwindows::toUnixSlash(currentDirectory);
    }
    else
    {
        FatalErrorInFunction
            << "Couldn't get the current working directory"
            << exit(FatalError);
    }

    return currentDirectory;
}


bool Foam::chDir(const fileName& dir)
{
    const bool success = ::SetCurrentDirectory(dir.c_str());
    return success;
}


bool Foam::mkDir(const fileName& pathName, const mode_t mode)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : pathName:" << pathName << " mode:" << mode
            << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (pathName.empty())
    {
        return false;
    }

    bool success = ::CreateDirectory(pathName.c_str(), nullptr);
    if (success)
    {
        chMod(pathName, mode);
    }
    else
    {
        const DWORD error = ::GetLastError();

        switch (error)
        {
            case ERROR_ALREADY_EXISTS:
            {
                success = true;
                break;
            }
            case ERROR_PATH_NOT_FOUND:
            {
                // Part of the path does not exist so try to create it
                const fileName& parentName = pathName.path();

                if (parentName.size() && mkDir(parentName, mode))
                {
                    success = mkDir(pathName, mode);
                }

                break;
            }
        }

        if (!success)
        {
            FatalErrorInFunction
              << "Couldn't create directory: " << pathName
              << " " << MSwindows::getLastError()
              << exit(FatalError);
        }
    }

    return success;
}


bool Foam::chMod(const fileName& name, const mode_t m)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    const int success = _chmod(name.c_str(), m);
    return success;
}


mode_t Foam::mode(const fileName& name, const bool followLink)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
    }

    fileStat fileStatus(name, followLink);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_mode;
    }
    else
    {
        return 0;
    }
}


Foam::fileName::Type Foam::type(const fileName& name, const bool followLink)
{
    //FIXME: 'followLink' can't be used with 'GetFileAttributes'.
    //Task for fixing this detail: https://github.com/blueCFD/Core/issues/60

    fileName::Type fileType = fileName::UNDEFINED;
    const DWORD attrs = ::GetFileAttributes(name.c_str());

    if (attrs != INVALID_FILE_ATTRIBUTES)
    {
        fileType = (attrs & FILE_ATTRIBUTE_DIRECTORY) ?
            fileName::DIRECTORY :
            fileName::FILE;
    }

    return fileType;
}


static
bool isGzFile(const std::string& name)
{
    const DWORD attrs = ::GetFileAttributes((name + ".gz").c_str());
    const bool success = (attrs != INVALID_FILE_ATTRIBUTES);

    return success;
}


bool Foam::exists
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    //FIXME: 'followLink' can't be used with 'GetFileAttributes'.
    //Task for fixing this detail: https://github.com/blueCFD/Core/issues/60

    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name
            << " checkGzip:" << checkGzip
            << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    const DWORD attrs = ::GetFileAttributes(name.c_str());
    const bool success = (attrs != INVALID_FILE_ATTRIBUTES) ||
                         (checkGzip && isGzFile(name));

    return success;
}


bool Foam::isDir(const fileName& name, const bool followLink)
{
    //FIXME: 'followLink' can't be used with 'GetFileAttributes'.
    //Task for fixing this detail: https://github.com/blueCFD/Core/issues/60

    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    const DWORD attrs = ::GetFileAttributes(name.c_str());
    bool success = (attrs != INVALID_FILE_ATTRIBUTES) &&
                   (attrs & FILE_ATTRIBUTE_DIRECTORY);

    return success;
}


bool Foam::isFile
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    //FIXME: 'followLink' can't be used with 'GetFileAttributes'.
    //Task for fixing this detail: https://github.com/blueCFD/Core/issues/60

   if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name
            << " checkGzip:" << checkGzip
            << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    const DWORD attrs = ::GetFileAttributes(name.c_str());
    const bool success =
    (
        (attrs != INVALID_FILE_ATTRIBUTES) &&
        !(attrs & FILE_ATTRIBUTE_DIRECTORY)
    )
     ||
    (
        checkGzip && isGzFile(name)
    );

    return success;
}


off_t Foam::fileSize(const fileName& name, const bool followLink)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    fileStat fileStatus(name, followLink);

    if (fileStatus.isValid())
    {
        return fileStatus.status().st_size;
    }
    else
    {
        return off_t(-1);
    }
}


time_t Foam::lastModified(const fileName& name, const bool followLink)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    fileStat fileStatus(name, followLink);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_mtime;
    }
    else
    {
        return 0;
    }
}


double Foam::highResLastModified(const fileName& name, const bool followLink)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    fileStat fileStatus(name, followLink);
    if (fileStatus.isValid())
    {
        return
            fileStatus.status().st_mtime
          + 1e-9*fileStatus.status().st_atim.tv_nsec;
    }
    else
    {
        return 0;
    }
}


Foam::fileNameList Foam::readDir
(
    const fileName& directory,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
)
{
    // Initial filename list size
    // also used as increment if initial size found to be insufficient
    const int maxNnames = 100;

    if (MSwindows::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : reading directory " << directory << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Setup empty string list MAXTVALUES long
    fileNameList dirEntries(maxNnames);

    // Temporary variables and counters
    label nEntries = 0;

    MSwindows::DirectoryIterator dirIt(directory);

    if (dirIt.isValid())
    {
        while (dirIt.hasNext())
        {
            const fileName & fName = dirIt.next();

            // ignore files begining with ., i.e. '.', '..' and '.*'
            if (fName.size() > 0 && fName[size_t(0)] != '.')
            {
                word fileNameExt = fName.ext();

                if
                (
                    (type == fileName::DIRECTORY)
                 ||
                    (
                        type == fileName::FILE
                        && fName[fName.size()-1] != '~'
                        && fileNameExt != "bak"
                        && fileNameExt != "BAK"
                        && fileNameExt != "old"
                        && fileNameExt != "save"
                    )
                )
                {
                    if ((directory/fName).type(followLink) == type)
                    {
                        if (nEntries >= dirEntries.size())
                        {
                            dirEntries.setSize(dirEntries.size() + maxNnames);
                        }

                        if (filtergz && fileNameExt == "gz")
                        {
                            dirEntries[nEntries++] = fName.lessExt();
                        }
                        else
                        {
                            dirEntries[nEntries++] = fName;
                        }
                    }
                }
            }
        }
    }
    else if (MSwindows::debug)
    {
        InfoInFunction
            << "cannot open directory " << directory << endl;
    }

    // Reset the length of the entries list
    dirEntries.setSize(nEntries);

    return dirEntries;
}


bool Foam::cp(const fileName& src, const fileName& dest, const bool followLink)
{
    if (MSwindows::debug)
    {
        Pout<< FUNCTION_NAME << " : src:" << src << " dest:" << dest << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Make sure source exists.
    if (!exists(src))
    {
        return false;
    }

    const fileName::Type srcType = src.type(followLink);

    fileName destFile(dest);

    // Check type of source file.
    if (srcType == fileName::FILE)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        // Open and check streams.
        // Use binary mode in case we read binary.
        // Causes windows reading to fail if we don't.
        std::ifstream srcStream(src.c_str(),
                                ios_base::in|ios_base::binary);
        if (!srcStream)
        {
            return false;
        }

        // Use binary mode in case we write binary.
        // Causes windows reading to fail if we don't.
        std::ofstream destStream(destFile.c_str(),
                                 ios_base::out|ios_base::binary);
        if (!destStream)
        {
            return false;
        }

        // Copy character data.
        char ch;
        while (srcStream.get(ch))
        {
            destStream.put(ch);
        }

        // Final check.
        if (!srcStream.eof() || !destStream)
        {
            return false;
        }
    }
    else if (srcType == fileName::LINK)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        ln(src, destFile);
    }
    else if (srcType == fileName::DIRECTORY)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.component(src.components().size() -1);
        }

        // Make sure the destination directory extists.
        if (!isDir(destFile) && !mkDir(destFile))
        {
            return false;
        }

        // Copy files
        fileNameList contents = readDir(src, fileName::FILE, false);
        forAll(contents, i)
        {
            if (MSwindows::debug)
            {
                InfoInFunction
                    << "Copying : " << src/contents[i]
                    << " to " << destFile/contents[i] << endl;
            }

            // File to file.
            cp(src/contents[i], destFile/contents[i]);
        }

        // Copy sub directories.
        fileNameList subdirs = readDir(src, fileName::DIRECTORY);
        forAll(subdirs, i)
        {
            if (MSwindows::debug)
            {
                InfoInFunction
                    << "Copying : " << src/subdirs[i]
                    << " to " << destFile << endl;
            }

            // Dir to Dir.
            cp(src/subdirs[i], destFile);
        }
    }

    return true;
}


bool Foam::ln(const fileName& src, const fileName& dest)
{
    // FIXME: Seems that prior to Vista softlinking was poorly supported.
    // Vista does a better job, but requires administrator privileges.
    // Task for this feature: https://github.com/blueCFD/Core/issues/63

    if (MSwindows::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : Create softlink from : " << src << " to " << dest << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (MSwindows::debug)
    {
        InfoInFunction
            << "MSwindows does not support ln - softlinking" << endl;
    }

    return false;
}


bool Foam::mv(const fileName& src, const fileName& dst, const bool followLink)
{
    if (MSwindows::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : Move : " << src << " to " << dst << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if
    (
        dst.type() == fileName::DIRECTORY
     && src.type(followLink) != fileName::DIRECTORY
    )
    {
        const fileName dstName(dst/src.name());

        return std::rename(src.c_str(), dstName.c_str()) == 0;
    }
    else
    {
        return std::rename(src.c_str(), dst.c_str()) == 0;
    }
}


bool Foam::mvBak(const fileName& src, const std::string& ext)
{
    if (MSwindows::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : moving : " << src << " to extension " << ext << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (exists(src, false))
    {
        const int maxIndex = 99;
        char index[3];

        for (int n = 0; n <= maxIndex; n++)
        {
            fileName dstName(src + "." + ext);
            if (n)
            {
                sprintf(index, "%02d", n);
                dstName += index;
            }

            // avoid overwriting existing files, except for the last
            // possible index where we have no choice
            if (!exists(dstName, false) || n == maxIndex)
            {
                return (0 == std::rename(src.c_str(), dstName.c_str()));
            }

        }
    }

    // fall-through: nothing to do
    return false;
}


// Remove a file returning true if successful otherwise false
bool Foam::rm(const fileName& file)
{
    if (MSwindows::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : Removing : " << file << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    bool success = (0 == std::remove(file.c_str()));

    // If deleting plain file name failed try with .gz
    if (!success)
    {
        const string fileGz = file + ".gz";
        success = (0 == std::remove(fileGz.c_str()));
    }

    return success;
}


// Remove a dirctory and it's contents
bool Foam::rmDir(const fileName& directory, const bool silent)
{
    if (MSwindows::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : removing directory " << directory << endl;
        if ((MSwindows::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    bool success = true;

    // Need to destroy DirectorIterator prior to
    // removing directory otherwise fails on Windows XP
    {
      MSwindows::DirectoryIterator dirIt(directory);

      while (success && dirIt.hasNext())
      {
          const fileName & fName = dirIt.next();

          if (fName != "." && fName != "..")
          {
              fileName path = directory/fName;

              if (path.type() == fileName::DIRECTORY)
              {
                  fileName const& const_path = path;
                  success = rmDir(const_path,true);

                  if (!success)
                  {
                      WarningInFunction
                        << "failed to remove directory " << fName
                        << " while removing directory " << directory
                        << endl;
                  }
              }
              else
              {
                  success = rm(path);

                  if (!success)
                  {
                      WarningInFunction
                        << "failed to remove file " << fName
                        << " while removing directory " << directory
                        << endl;
                  }
              }
          }
      }
    }

    if (success)
    {
        success = ::RemoveDirectory(directory.c_str());

        if (!success)
        {
            WarningInFunction
                << "failed to remove directory " << directory << endl;
        }
    }

    return success;
}


unsigned int Foam::sleep(const unsigned int s)
{
    const DWORD milliseconds = s * 1000;

    ::Sleep(milliseconds);

    return 0;
}

int Foam::nanosleep(const scalar& s)
{
    // MINGW does not have the nanosleep function
    // use a modified sleep function
    const DWORD nanoseconds = s * 1;

    ::Sleep(nanoseconds);

    return 0;
}


void Foam::fdClose(const int fd)
{
    const int result = ::_close(fd);

    if (0 != result)
    {
        FatalErrorInFunction
            << "close error on " << fd << endl
            << abort(FatalError);
    }
}


bool Foam::ping
(
    const std::string& destName,
    const label destPort,
    const label timeOut
)
{
    // FIXME: Task for this: https://github.com/blueCFD/Core/issues/64

    if (MSwindows::debug)
    {
        InfoInFunction
            << "MSwindows does not support ping" << endl;
    }

    return false;
}


bool Foam::ping(const std::string&  hostname, const label timeOut)
{
    return ping(hostname, 222, timeOut) || ping(hostname, 22, timeOut);
}


int Foam::system(const std::string& command)
{
    return std::system(command.c_str());
}

int Foam::system(const Foam::UList<Foam::string>& command)
{
    return std::system(command[0].c_str());
}

// Explicitly track loaded libraries, rather than use
// EnumerateLoadedModules64 and have to link against Dbghelp.dll
// Details at http://msdn.microsoft.com/en-us/library/ms679316(v=vs.85).aspx
typedef std::map<void*, std::string> OfLoadedLibs;

static
OfLoadedLibs &
getLoadedLibs()
{
  static OfLoadedLibs loadedLibs;
  return loadedLibs;
}


void* Foam::dlOpen(const fileName& libName, const bool check)
{
    //Lets check if this is a list of libraries to be loaded
    //NOTE: should only be used for "force loading libraries"
    if (libName.find_first_of(',')!=Foam::string::npos)
    {
        void *moduleh=nullptr;
        string libsToLoad=libName;
        libsToLoad.removeTrailing(' '); //removes spaces from both ends
        libsToLoad.removeRepeated(',');
        libsToLoad += ',';

        if (MSwindows::debug)
        {
            InfoInFunction
                << "Libraries to be loaded: " <<  libsToLoad << endl;
        }

        //generate the word list
        size_t stposstr=0, found=libsToLoad.find_first_of(',');
        while (found!=string::npos)
        {
            string libToLoad = libsToLoad.substr(stposstr,found-stposstr);
            moduleh = dlOpen(libToLoad);
            stposstr=found+1; found=libsToLoad.find_first_of(',',stposstr);
        }

        return moduleh;
    }
    else
    {
        if (MSwindows::debug)
        {
            InfoInFunction
                << " : LoadLibrary of " << libName << endl;
        }

        const char* dllExt = ".dll";

        // Assume libName is of the form, lib<name>.so
        Foam::string winLibName(libName);
        winLibName.replace(".so", dllExt);
        void* libHandle = ::LoadLibrary(winLibName.c_str());

        if (nullptr == libHandle)
        {
            // Assumes libName = name
            winLibName = "lib";
            winLibName += libName;
            winLibName += dllExt;

            libHandle = ::LoadLibrary(winLibName.c_str());
        }

        if (nullptr == libHandle && check)
        {
            WarningInFunction
                << "Failed to load library with name "
                << libName << ": "
                << MSwindows::getLastError()
                << endl;
        }
        else
        {
            getLoadedLibs()[libHandle] = libName;
        }

        if (MSwindows::debug)
        {
            InfoInFunction
                << "Library " <<  libName << " loaded "
                << (libHandle != nullptr ? "with success!" : "without success.")
                << endl;
        }

        return libHandle;
    }
}


bool Foam::dlClose(void* libHandle)
{
    if (MSwindows::debug)
    {
        InfoInFunction
            << "FreeLibrary of handle " << libHandle << endl;
    }

    const bool success =
      ::FreeLibrary(static_cast<HMODULE>(libHandle));

    if (!success)
    {
        WarningInFunction
            << "FreeLibrary failed. "
            << MSwindows::getLastError()
            << endl;
    }
    else
    {
        getLoadedLibs().erase(libHandle);
    }

    return success;
}


void* Foam::dlSym(void* handle, const std::string& symbol)
{
    if (MSwindows::debug)
    {
        InfoInFunction
            << "dlsym of " << symbol << endl;
    }

    // get address of symbol
    void* fun = (void *)
    (
        ::GetProcAddress(static_cast<HMODULE>(handle), symbol.c_str())
    );

    if (fun == nullptr)
    {
        WarningInFunction
            << "Cannot lookup symbol " << symbol << " : "
            << MSwindows::getLastError()
            << endl;
    }

    return fun;
}


bool Foam::dlSymFound(void* handle, const std::string& symbol)
{
    if (handle && !symbol.empty())
    {
        if (MSwindows::debug)
        {
            InfoInFunction
                << "dlSym of " << symbol << endl;
        }

        // symbol can be found if there was no error
        return dlSym(handle, symbol.c_str()) != nullptr;
    }
    else
    {
        return false;
    }
}


Foam::fileNameList Foam::dlLoaded()
{
    fileNameList libs;
    OfLoadedLibs & loadedLibs = getLoadedLibs();

    for (OfLoadedLibs::const_iterator it = loadedLibs.begin();
         it != loadedLibs.end();
         ++it)
    {
        libs.append(it->second);
    }

    if (MSwindows::debug)
    {
        InfoInFunction
            << "determined loaded libraries :" << libs.size() << endl;
    }
    return libs;
}

//- Random functions

//It's easier to include it here...
#include "random.c"

void Foam::osRandomSeed(const label seed)
{
    srandom((unsigned int)(seed));
}

Foam::label Foam::osRandomInteger()
{
    return random();
}

Foam::scalar Foam::osRandomDouble()
{
    return scalar(random()) / scalar(INT_MAX);
}

Foam::string Foam::toUnixPath(const string & path)
{
    string unixPath(path);
    MSwindows::toUnixSlash(unixPath);
    MSwindows::removeQuotes(unixPath);

    return unixPath;
}


//- Thread handling: Using std::thread and std::mutex

static Foam::DynamicList<Foam::autoPtr<std::thread>> threads_;
static Foam::DynamicList<Foam::autoPtr<std::mutex>> mutexes_;


Foam::label Foam::allocateThread()
{
    forAll(threads_, i)
    {
        if (!threads_[i].valid())
        {
            if (MSwindows::debug)
            {
                Pout<< "allocateThread : reusing index:" << i << endl;
            }
            // Reuse entry
            threads_[i].reset(new std::thread());
            return i;
        }
    }

    label index = threads_.size();
    if (MSwindows::debug)
    {
        Pout<< "allocateThread : new index:" << index << endl;
    }
    threads_.append(autoPtr<std::thread>(new std::thread()));

    return index;
}


void Foam::createThread
(
    const label index,
    void *(*start_routine) (void *),
    void *arg
)
{
    if (MSwindows::debug)
    {
        Pout<< "createThread : index:" << index << endl;
    }

    try
    {
        threads_[index].reset(new std::thread(start_routine, arg));
    }
    catch(...)
    {
        FatalErrorInFunction
            << "Failed starting thread " << index << exit(FatalError);
    }
}


void Foam::joinThread(const label index)
{
    if (MSwindows::debug)
    {
        Pout<< "joinThread : join:" << index << endl;
    }

    if (!threads_[index]().joinable())
    {
        FatalErrorInFunction << "Failed freeing thread " << index
            << exit(FatalError);
    }
    else
    {
        try
        {
            threads_[index]().join();
        }
        catch(...)
        {
            FatalErrorInFunction << "Failed freeing thread " << index
                << exit(FatalError);
        }
    }
}


void Foam::freeThread(const label index)
{
    if (MSwindows::debug)
    {
        Pout<< "freeThread : index:" << index << endl;
    }
    threads_[index].clear();
}

Foam::label Foam::allocateMutex()
{
    forAll(mutexes_, i)
    {
        if (!mutexes_[i].valid())
        {
            if (MSwindows::debug)
            {
                Pout<< "allocateMutex : reusing index:" << i << endl;
            }
            // Reuse entry
            mutexes_[i].reset(new std::mutex());
            return i;
        }
    }

    label index = mutexes_.size();

    if (MSwindows::debug)
    {
        Pout<< "allocateMutex : new index:" << index << endl;
    }
    mutexes_.append(autoPtr<std::mutex>(new std::mutex()));
    return index;
}


void Foam::lockMutex(const label index)
{
    if (MSwindows::debug)
    {
        Pout<< "lockMutex : index:" << index << endl;
    }

    try
    {
        mutexes_[index]().lock();
    }
    catch(...)
    {
        FatalErrorInFunction << "Failed locking mutex " << index
            << exit(FatalError);
    }
}


void Foam::unlockMutex(const label index)
{
    if (MSwindows::debug)
    {
        Pout<< "unlockMutex : index:" << index << endl;
    }

    try
    {
        mutexes_[index]().unlock();
    }
    catch(...)
    {
        FatalErrorInFunction << "Failed unlocking mutex " << index
            << exit(FatalError);
    }
}


void Foam::freeMutex(const label index)
{
    if (MSwindows::debug)
    {
        Pout<< "freeMutex : index:" << index << endl;
    }
    mutexes_[index].clear();
}

// ************************************************************************* //
