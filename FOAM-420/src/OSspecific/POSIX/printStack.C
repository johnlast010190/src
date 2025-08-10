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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "db/IOstreams/StringStreams/OStringStream.H"
#include "include/OSspecific.H"
#include "db/IOstreams/Fstreams/IFstream.H"

#include <inttypes.h>
#include <cxxabi.h>
#include <execinfo.h>
#include <dlfcn.h>
#include <unistd.h>
#ifdef USE_BFD
    // Some versions of binutils require these variables to be set.  This works
    // around the following bug:
    // bfd.h(35): error: #error directive: config.h must be included before this header
    // This bug seen with binutils version 2.28, but binutils v2.37 seems to
    // work fine.
    #ifndef PACKAGE
        #define PACKAGE
        #define UNDEF_PACKAGE
    #endif
    #ifndef PACKAGE_VERSION
        #define PACKAGE_VERSION
        #define UNDEF_PACKAGE_VERSION
    #endif
    #include <bfd.h>
    #ifdef UNDEF_PACKAGE
        #undef PACKAGE
        #undef UNDEF_PACKAGE
    #endif
    #ifdef UNDEF_PACKAGE_VERSION
        #undef PACKAGE_VERSION
        #undef UNDEF_PACKAGE_VERSION
    #endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string pOpen(const string &cmd, label line=0)
{
    string res = "\n";

    FILE *cmdPipe = popen(cmd.c_str(), "r");

    if (cmdPipe)
    {
        char *buf = nullptr;

        // Read line number of lines
        for (label cnt = 0; cnt <= line; cnt++)
        {
            size_t linecap = 0;
            ssize_t linelen;
            linelen = getline(&buf, &linecap, cmdPipe);

            if (linelen < 0)
            {
                break;
            }

            if (cnt == line)
            {
                res = string(buf);
                break;
            }
        }

        if (buf != nullptr)
        {
            free(buf);
        }

        pclose(cmdPipe);
    }

    return res.substr(0, res.size() - 1);
}


inline word addressToWord(const uintptr_t addr)
{
    OStringStream nStream;
    nStream << "0x" << hex << addr;
    return nStream.str();
}


bool addrToLine
(
    const fileName& filename,
    const intptr_t relativeAddress,
    string& outputLine
)
{
    if (env("FOAM_USE_ADDR2LINE"))
    {
        // Force the old behaviour - not recommended as forking can cause
        // trouble with certain OpenMPI drivers
        word myAddress = addressToWord(relativeAddress);
        outputLine = pOpen
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + myAddress,
            1
        );
        if (outputLine == "" || outputLine == "??:?" || outputLine == "??:0")
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        #ifdef USE_BFD

            #if not defined(bfd_get_section_name)
            #define BFD_2_34
            #endif

            bfd_init();

            bfd_vma addr = relativeAddress;

            bfd* abfd = bfd_openr(filename.c_str(), NULL);
            if (!abfd)
            {
                return false;
            }
            char** matching;
            if
            (
                bfd_check_format(abfd, bfd_archive)
             || !bfd_check_format_matches(abfd, bfd_object, &matching)
             || !(bfd_get_file_flags(abfd) & HAS_SYMS)
            )
            {
                bfd_close(abfd);
                return false;
            }

            asymbol** syms;
            unsigned int size;
            long symcount =
                bfd_read_minisymbols
                (
                    abfd, false, reinterpret_cast<void**>(&syms), &size
                );
            if (symcount == 0)
            {
                symcount =
                    bfd_read_minisymbols
                    (
                        abfd, true, reinterpret_cast<void**>(&syms), &size
                    );
            }
            if (symcount < 0)
            {
                bfd_close(abfd);
                return false;
            }

            // Preferred method is to use bfd_map_over_sections, but this loop
            // is more compact
            asection *sect;
            bool found = false;
            for (sect = abfd->sections; sect != nullptr; sect = sect->next)
            {
                #ifdef BFD_2_34
                if ((bfd_section_flags(sect) & SEC_ALLOC) == 0)
                #else
                if ((bfd_get_section_flags(abfd, sect) & SEC_ALLOC) == 0)
                #endif
                {
                    bfd_close(abfd);
                    free(syms);
                    return false;
                }

                bfd_vma vma;
                #ifdef BFD_2_34
                vma = bfd_section_vma(sect);
                #else
                vma = bfd_get_section_vma(abfd, sect);
                #endif
                if (addr < vma)
                {
                    bfd_close(abfd);
                    free(syms);
                    return false;
                }

                bfd_size_type size;
                #ifdef BFD_2_34
                size = bfd_section_size(sect);
                #else
                size = bfd_get_section_size(sect);
                #endif
                if (addr >= vma + size)
                {
                    continue;
                }

                const char *sourceFile = nullptr;
                const char *method = nullptr;
                unsigned int line = 0;
                found =
                    bfd_find_nearest_line
                    (
                        abfd,
                        sect,
                        syms,
                        addr - vma,
                        &sourceFile,
                        &method,
                        &line
                    );

                if (found)
                {
                    std::string sourceFileName;
                    if (sourceFile == nullptr)
                    {
                        found = false;
                    }
                    else
                    {
                        sourceFileName = sourceFile;

                        OStringStream oss;
                        oss << sourceFile << ":";
                        if (line > 0)
                        {
                            oss << line;
                        }
                        else
                        {
                            oss << "?";
                        }
                        outputLine = oss.str();
                    }

                    break;
                }
            }
            free(syms);
            bfd_close(abfd);

            return found;

        #else

            return false;

        #endif
    }
}


void printSourceFileAndLine
(
    Ostream& os,
    const fileName& filename,
    Dl_info *info,
    void *addr
)
{
    uintptr_t relativeAddress = uintptr_t(addr);
    if (filename.ext() == "so")
    {
        // Convert address into offset into dynamic library
        uintptr_t offset = uintptr_t(info->dli_fbase);
        relativeAddress -= offset;
    }

    if (filename[0] == '/')
    {
        string line;
        bool found = addrToLine(filename, relativeAddress, line);

        if (!found)
        {
            os  << " in " << filename;
        }
        else
        {
            string cwdLine(line.replaceAll(cwd() + '/', ""));
            string homeLine(cwdLine.replaceAll(home(), '~'));

            os  << " at " << homeLine.c_str();
        }
    }
}


fileName absolutePath(const char* fn)
{
    fileName fname(fn);

    if (fname[0] != '/' && fname[0] != '~')
    {
        string tmp;
        if (env("FOAM_USE_ADDR2LINE"))
        {
            // Force the old behaviour - not recommended as forking can cause
            // trouble with certain OpenMPI drivers
            tmp = pOpen("which " + fname);
        }
        else
        {
            string path = getEnv("PATH");
            while (true)
            {
                size_t pos = path.find(":");
                string fullPath = path.substr(0, pos);
                if (fullPath.size())
                {
                    if (fullPath[fullPath.size()-1] != '/')
                    {
                        fullPath += '/';
                    }
                    fullPath += fname;
                    if (access(fullPath.c_str(), X_OK) == 0)
                    {
                        tmp = fullPath;
                        break;
                    }
                }
                if (pos == string::npos)
                {
                    break;
                }
                path = path.substr(pos+1);
            }
        }

        if (tmp[0] == '/' || tmp[0] == '~')
        {
            fname = tmp;
        }
    }

    return fname;
}


word demangleSymbol(const char* sn)
{
    word res;
    int st;
    char* cxx_sname = abi::__cxa_demangle
    (
        sn,
        nullptr,
        0,
        &st
    );

    if (st == 0 && cxx_sname)
    {
        res = word(cxx_sname);
        free(cxx_sname);
    }
    else
    {
        res = word(sn);
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::error::safePrintStack(std::ostream& os)
{
    // Get raw stack symbols
    void *array[100];
    size_t size = backtrace(array, 100);
    char **strings = backtrace_symbols(array, size);

    // See if they contain function between () e.g. "(__libc_start_main+0xd0)"
    // and see if cplus_demangle can make sense of part before +
    for (size_t i = 0; i < size; i++)
    {
        string msg(strings[i]);
        fileName programFile;
        word address;

        os  << '#' << label(i) << '\t' << msg << std::endl;
    }
}


void Foam::error::printStack(Ostream& os)
{
    // Get raw stack symbols
    const size_t CALLSTACK_SIZE = 128;

    void *callstack[CALLSTACK_SIZE];
    size_t size = backtrace(callstack, CALLSTACK_SIZE);

    Dl_info *info = new Dl_info;

    fileName fname = "???";
    word address;

    for (size_t i=0; i<size; i++)
    {
        int st = dladdr(callstack[i], info);

        os << '#' << label(i) << "  ";
        if (st != 0 && info->dli_fname != nullptr && info->dli_fname[0] != '\0')
        {
            fname = absolutePath(info->dli_fname);

            os <<
            (
                (info->dli_sname != nullptr)
              ? demangleSymbol(info->dli_sname)
              : "?"
            );
        }
        else
        {
            os << "?";
        }

        printSourceFileAndLine(os, fname, info, callstack[i]);
        os << nl;
    }

    delete info;
}


// ************************************************************************* //
