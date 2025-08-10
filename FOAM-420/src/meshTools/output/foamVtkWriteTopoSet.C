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
    (c) 2018 OpenCFD Ltd.


\*---------------------------------------------------------------------------*/

#include "output/foamVtkWriteTopoSet.H"
#include "meshes/polyMesh/polyMesh.H"
#include "sets/topoSets/topoSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/pointSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::vtk::writeTopoSet
(
    const polyMesh& mesh,
    const topoSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
{
    if (isA<pointSet>(set))
    {
        return vtk::writePointSet
        (
            mesh,
            dynamicCast<const pointSet&>(set),
            opts,
            file,
            parallel
        );
    }
    else if (isA<faceSet>(set))
    {
        return vtk::writeFaceSet
        (
            mesh,
            dynamicCast<const faceSet&>(set),
            opts,
            file,
            parallel
        );
    }
    else if (isA<cellSet>(set))
    {
        return vtk::writeCellSetFaces
        (
            mesh,
            dynamicCast<const cellSet&>(set),
            opts,
            file,
            parallel
        );
    }

    WarningInFunction
        << "No VTK writer for '" << set.type() << "' topoSet" << nl << endl;

    return false;
}


// ************************************************************************* //
