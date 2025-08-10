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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.


Description

\*---------------------------------------------------------------------------*/

#include "faMesh/faPatches/constraint/wedge/wedgeFaPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "faMesh/faBoundaryMesh/faBoundaryMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.H"
#include "meshes/polyMesh/polyMesh.H"
#include "faMesh/faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(wedgeFaPatch, 0);
addToRunTimeSelectionTable(faPatch, wedgeFaPatch, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void wedgeFaPatch::findAxisPoint() const
{
    // Find axis point

    labelList ptLabels = pointLabels();

    labelListList ptEdges = pointEdges();

    const vectorField& points = boundaryMesh().mesh().points();

    const scalarField& magL = magEdgeLengths();

    forAll(ptEdges, pointI)
    {
        if (ptEdges[pointI].size() == 1)
        {
            scalar r = mag((I-axis()*axis())&points[ptLabels[pointI]]);

            if (r < magL[ptEdges[pointI][0]])
            {
                axisPoint_ = ptLabels[pointI];
                break;
            }
        }
    }

    axisPointChecked_ = true;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

//- Construct from polyPatch
wedgeFaPatch::wedgeFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    faPatch(name, dict, index, bm),
    wedgePolyPatchPtr_(nullptr),
    axisPoint_(-1),
    axisPointChecked_(false)
{
    if (ngbPolyPatchIndex() == -1)
    {
        FatalErrorInFunction
            << "Neighbour polyPatch index is not specified for faPatch "
            << this->name()
            << exit(FatalError);
    }

    if
    (
        bm.mesh()().boundaryMesh()[ngbPolyPatchIndex()].type()
     == wedgePolyPatch::typeName
    )
    {
        const wedgePolyPatch& wedge =
            refCast<const wedgePolyPatch>
            (
                bm.mesh()().boundaryMesh()[ngbPolyPatchIndex()]
            );

        wedgePolyPatchPtr_ = &wedge;
    }
    else
    {
        FatalErrorInFunction
            << "Neighbour polyPatch is not of type "
            << wedgePolyPatch::typeName
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
