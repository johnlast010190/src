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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "patchRotatingVelocityForce.H"
#include "meshes/polyMesh/polyPatches/constraint/processorCyclic/processorCyclicPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(patchRotatingVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        patchRotatingVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::patchRotatingVelocityForce::patchRotatingVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    rotatingVelocityForce(sourceName, modelType, dict, obr),
    patch_(coeffs_.lookup("patch")),
    patchi_(mesh_.boundaryMesh().findPatchID(patch_))
{
    if (patchi_ < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << patch_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::patchRotatingVelocityForce::omegaAve
(
    const volVectorField& U
) const
{
    const fvPatch& cPatch(mesh_.boundary()[patchi_]);
    const vectorField& Up(U.boundaryField()[patchi_]);
    vectorField axR(cPatch.Cf() - origin_);
    axR -= axis_*(axR & axis_);
    axR = (axis_ ^ axR)();

    vector2D sumAomegaSumA
    (
        sum(((Up & axR) / magSqr(axR))*cPatch.magSf()),
        sum(cPatch.magSf())
    );


    // If the rotating velocity force is applied to a cyclic patch
    // for parallel runs include contributions from processorCyclic patches
    // generated from the decomposition of the cyclic patch
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun() && isA<cyclicPolyPatch>(patches[patchi_]))
    {
        labelList processorCyclicPatches
        (
            processorCyclicPolyPatch::patchIDs(patch_, patches)
        );

        forAll(processorCyclicPatches, pcpi)
        {
            const label patchi = processorCyclicPatches[pcpi];


            const fvPatch& pcPatch(mesh_.boundary()[patchi]);
            const vectorField& Upc(U.boundaryField()[patchi]);
            axR = (pcPatch.Cf() - origin_);
            axR -= axis_*(axR & axis_);
            axR = (axis_ ^ axR);

            sumAomegaSumA.x()
                += sum(((Upc & axR) / magSqr(axR))*pcPatch.magSf());

            sumAomegaSumA.y() += sum(pcPatch.magSf());
        }
    }

    mesh_.reduce(sumAomegaSumA, sumOp<vector2D>());

    return sumAomegaSumA.x()/sumAomegaSumA.y();
}


// ************************************************************************* //
