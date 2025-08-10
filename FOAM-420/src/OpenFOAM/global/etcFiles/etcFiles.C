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
    (c) 2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "global/etcFiles/etcFiles.H"
#include "include/OSspecific.H"
#include "global/foamVersion.H"

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

//
// These could be exposed too (if required), but are fairly special purpose.
//
//! \cond fileScope
//
// Assign 'queried' parameter to the user resource directory.
// Return true if this directory exists.
//
// Corresponds to foamEtcFile -mode=u
// Looks for
//   - ~/.FOAMcore
static inline bool userResourceDir(Foam::fileName& queried)
{
    queried = Foam::home()/FOAM_USER_RESOURCE_DIRNAME;
    return Foam::isDir(queried);
}


// Assign 'queried' parameter to the group resource directory.
// Return true if this directory exists.
//
// Corresponds to foamEtcFile -mode=g
// Looks for
//   - $FOAM_CUSTOM_DATA_DIR
//   - FOAM_PROJECT_DIR/../site
static inline bool groupResourceDir(Foam::fileName& queried)
{
    queried = Foam::getEnv("FOAM_CUSTOM_DATA_DIR");
    if (queried.size())
    {
        return Foam::isDir(queried);
    }

    // Fallback (when FOAM_CUSTOM_DATA_DIR is unset)
    queried = Foam::getEnv("FOAM_PROJECT_DIR")/".."/"site";
    return (queried.size() > 5 && Foam::isDir(queried));
}


// Assign 'queried' parameter to the OpenFOAM etc/ resource directory.
// Return true if it exists.
//
// Corresponds to foamEtcFile -mode=o
// Looks for
//   - $FOAM_PROJECT_DIR/etc
static inline bool projectResourceDir(Foam::fileName& queried)
{
    queried = Foam::getEnv("FOAM_PROJECT_DIR")/"etc";
    return (queried.size() > 4 && Foam::isDir(queried));
}

//! \endcond


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileNameList Foam::findEtcDirs
(
    const fileName& name,
    const bool findFirst
)
{
    fileNameList results;

    do
    {
        fileName dir, candidate;

        // User resource directories
        if (userResourceDir(dir))
        {
            candidate = dir/FOAMversion/name;
            if (isDir(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }

            candidate = dir/name;
            if (isDir(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }
        }

        // Group resource directories
        if (groupResourceDir(dir))
        {
            candidate = dir/FOAMversion/name;
            if (isDir(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }

            candidate = dir/name;
            if (isDir(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }
        }

        // Other (project) resource directory
        if (projectResourceDir(dir))
        {
            candidate = dir/name;
            if (isDir(dir) && isDir(candidate))
            {
                results.append(candidate);
            }
        }
    }
    while (false);  // Run exactly once

    return results;
}


Foam::fileNameList Foam::findEtcFiles
(
    const fileName& name,
    const bool mandatory,
    const bool findFirst
)
{
    fileNameList results;

    do
    {
        fileName dir, candidate;

        // User resource directories
        if (userResourceDir(dir))
        {
            candidate = dir/FOAMversion/name;
            if (isFile(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }

            candidate = dir/name;
            if (isFile(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }
        }

        // Group resource directories
        if (groupResourceDir(dir))
        {
            candidate = dir/FOAMversion/name;
            if (isFile(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }

            candidate = dir/name;
            if (isFile(candidate))
            {
                results.append(candidate);
                if (findFirst)
                {
                    break;
                }
            }
        }

        // Other (project) resource directory
        if (projectResourceDir(dir))
        {
            candidate = dir/name;
            if (isDir(dir) && isFile(candidate))
            {
                results.append(candidate);
            }
        }
    }
    while (false);  // Run exactly once

    // No name?  It cannot be a file!
    if (name.empty())
    {
        results.clear();
    }

    if (mandatory && results.empty())
    {
        // Abort if file is mandatory but not found
        std::cerr
            << "--> FOAM FATAL ERROR in Foam::findEtcFiles()"
               " :  could not find mandatory file\n    '"
            << name.c_str() << "'\n\n" << std::endl;
        ::exit(1);
    }

    return results;
}


Foam::fileName Foam::findEtcFile(const fileName& name, const bool mandatory)
{
    fileNameList results(findEtcFiles(name, mandatory, true));

    if (results.size())
    {
        return results[0];
    }

    // Return null-constructed fileName rather than fileName::null
    return fileName();
}


// ************************************************************************* //
