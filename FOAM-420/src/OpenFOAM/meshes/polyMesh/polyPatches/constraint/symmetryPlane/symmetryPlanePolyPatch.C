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
    (c) 2013-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyPatches/constraint/symmetryPlane/symmetryPlanePolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyPatches/constraint/symmetry/symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symmetryPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, dictionary);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::symmetryPlanePolyPatch::calcGeometry(PstreamBuffers&)
{
    if (n_ == vector::rootMax)
    {
        if (returnReduce(size(), sumOp<label>()))
        {
            const vectorField& nf(faceNormals());
            n_ = gAverage(nf);

            if (debug)
            {
                Info<< "Patch " << name() << " calculated average normal "
                    << n_ << endl;
            }


            // Check the symmetry plane is planar
            forAll(nf, facei)
            {
                if (magSqr(n_ - nf[facei]) > SMALL)
                {
                    WarningInFunction
                        << "Symmetry plane '" << name() << "' is not planar."
                        << endl
                        << "At local face at "
                        << primitivePatch::faceCentres()[facei]
                        << " the normal " << nf[facei]
                        << " differs from the average normal " << n_
                        << " by " << magSqr(n_ - nf[facei]) << endl
                        << "Either split the patch into planar parts"
                        << " or use the " << symmetryPolyPatch::typeName
                        << " patch type"
                        << endl;//exit(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    directPolyPatch(name, size, start, index, bm, patchType),
    n_(vector::rootMax)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    directPolyPatch(name, dict, index, bm, patchType),
    n_(vector::rootMax)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    directPolyPatch(pp, bm),
    n_(pp.n_)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, newSize, newStart),
    n_(pp.n_)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, mapAddressing, newStart),
    n_(pp.n_)
{}


// ************************************************************************* //
