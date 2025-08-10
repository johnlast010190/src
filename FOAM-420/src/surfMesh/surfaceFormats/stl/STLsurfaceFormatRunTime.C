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
    (c) 2011 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFormats/stl/STLsurfaceFormat.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

// read MeshedSurface (ascii)
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STLsurfaceFormat,
    face,
    fileExtension,
    stl
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STLsurfaceFormat,
    triFace,
    fileExtension,
    stl
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STLsurfaceFormat,
    labelledTri,
    fileExtension,
    stl
);

// read MeshedSurface (binary)
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STLsurfaceFormat,
    face,
    fileExtension,
    stlb
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STLsurfaceFormat,
    triFace,
    fileExtension,
    stlb
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STLsurfaceFormat,
    labelledTri,
    fileExtension,
    stlb
);


// write MeshedSurfaceProxy (ascii)
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STLsurfaceFormat,
    face,
    write,
    fileExtension,
    stl
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STLsurfaceFormat,
    triFace,
    write,
    fileExtension,
    stl
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STLsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    stl
);

// write MeshedSurfaceProxy (binary)
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STLsurfaceFormat,
    face,
    write,
    fileExtension,
    stlb
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STLsurfaceFormat,
    triFace,
    write,
    fileExtension,
    stlb
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STLsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    stlb
);

// write UnsortedMeshedSurface (ascii)
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    STLsurfaceFormat,
    face,
    write,
    fileExtension,
    stl
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    STLsurfaceFormat,
    triFace,
    write,
    fileExtension,
    stl
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    STLsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    stl
);

// write UnsortedMeshedSurface (binary)
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    STLsurfaceFormat,
    face,
    write,
    fileExtension,
    stlb
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    STLsurfaceFormat,
    triFace,
    write,
    fileExtension,
    stlb
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    UnsortedMeshedSurface,
    STLsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    stlb
);

}
}

// ************************************************************************* //
