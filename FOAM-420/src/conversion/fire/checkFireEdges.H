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
    (c) 2016 OpenCFD Ltd.

Description
    Checks the mesh for edge connectivity as expected by the AVL/FIRE
    on-the-fly calculations.
    Errors flagged here are not necessarily topological errors at all.

SourceFiles
    checkFireEdges.C

\*---------------------------------------------------------------------------*/

#ifndef checkFireEdges_H
#define checkFireEdges_H

#include "meshes/meshShapes/face/faceList.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;

//- check edge connectivity
label checkFireEdges
(
    const faceList& faces,
    const labelListList& pointFaces,
    const UList<point>& points = UList<point>::null()
);


//- check edge connectivity with pointFaces mapping calculated automatically
label checkFireEdges
(
    const faceList&,
    const UList<point>& points = UList<point>::null()
);


//- check edge connectivity
label checkFireEdges(const polyMesh& mesh);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
