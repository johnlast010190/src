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

InNamespace
    Foam

Description
    Functions to search 'etc' directories for configuration files etc.

SourceFiles
    etcFiles.C

\*---------------------------------------------------------------------------*/

#ifndef etcFiles_H
#define etcFiles_H

#include "primitives/strings/lists/fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Search for directories from user/group/other directories.
//  The search scheme allows for version-specific and
//  version-independent files using the following hierarchy:
//  - \b user settings:
//    - ~/.FOAMcore/\<VERSION\>
//    - ~/.FOAMcore/
//  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is set):
//    - $FOAM_CUSTOM_DATA_DIR/\<VERSION\>
//    - $FOAM_CUSTOM_DATA_DIR
//  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is not set):
//    - $FOAM_PROJECT_DIR/../site/\<VERSION\>
//    - $FOAM_PROJECT_DIR/../site/
//  - \b other (project) settings:
//    - $FOAM_PROJECT_DIR/etc/
//
//  \return The list of full paths of all the matching directories or
//  an empty list if the name cannot be found.
//  Optionally stop search after the first directory has been found.
fileNameList findEtcDirs
(
    const fileName& name = fileName::null,
    const bool findFirst = false
);

//- Search for files from user/group/other directories.
//  The search scheme allows for version-specific and
//  version-independent files using the following hierarchy:
//  - \b user settings:
//    - ~/.FOAMcore/\<VERSION\>
//    - ~/.FOAMcore/
//  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is set):
//    - $FOAM_CUSTOM_DATA_DIR/\<VERSION\>
//    - $FOAM_CUSTOM_DATA_DIR
//  - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is not set):
//    - $FOAM_PROJECT_DIR/../site/\<VERSION\>
//    - $FOAM_PROJECT_DIR/../site/
//  - \b other (project) settings:
//    - $FOAM_PROJECT_DIR/etc/
//
//  \return The list of full paths of all the matching files or
//  an empty list if the name cannot be found.
//  Optionally abort if the file cannot be found.
//  Optionally stop search after the first file has been found.
fileNameList findEtcFiles
(
    const fileName& name,
    const bool mandatory = false,
    const bool findFirst = false
);

//- Search for a single file using findEtcFiles.
//  \return The full path name of the first file found in the
//  search hierarchy or an empty fileName if the name cannot be found.
//  Optionally abort if the file cannot be found but is mandatory.
fileName findEtcFile(const fileName& name, const bool mandatory=false);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
