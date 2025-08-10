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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/adjustU/adjustU.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchFields.H"
#include "fields/fvPatchFields/derived/windProfileDirectionVelocity/windProfileDirectionVelocityFvPatchVectorField.H"
#include "fields/fvPatchFields/derived/outflowVelocity/outflowVelocityFvPatchVectorField.H"
#include "fields/fvPatchFields/derived/interpolatedInletOutlet/interpolatedInletOutletFvPatchFields.H"
#include "fields/fvPatchFields/basic/blended/blendedFvPatchFields.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::adjustU
(
    const surfaceScalarField& phi,
    const word rhoName,
    volVectorField& U
)
{
    bool outflowExist = false;
    forAll(U.boundaryField(), pI)
    {
        const fvPatchVectorField& Up = U.boundaryField()[pI];
        if
        (
            isA<outflowVelocityFvPatchVectorField>(Up)
        )
        {
            outflowExist = true;
        }
    }

    if (outflowExist)
    {
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;
        scalar areaAdjustable = 0.0;
        scalar rhoAdjustable = 0.0;

        bool compressible((phi.dimensions() == dimDensity*dimVelocity*dimArea));

        const surfaceScalarField::Boundary& bphi =
            phi.boundaryField();
        const fvMesh& mesh = phi.mesh();

        forAll(bphi, patchi)
        {
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];
            const scalarField& patchAreas = mesh.boundary()[patchi].magSf();

            if (!phip.coupled())
            {
                if
                (
                    Up.fixesValue()
                 && !isA<inletOutletFvPatchVectorField>(Up)
                 && !isA<interpolatedInletOutletFvPatchVectorField>(Up)
                 && !isA<windProfileDirectionVelocityFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }

                }
                else
                {
                    if
                    (
                        isA<outflowVelocityFvPatchVectorField>(Up)
                    )
                    {
                        forAll(phip, i)
                        {
                            adjustableMassOut += phip[i];
                        }
                        areaAdjustable += sum(patchAreas);
                        if (compressible)
                        {

                            const volScalarField& rho =
                                U.db().lookupObject<volScalarField>
                                (
                                    rhoName
                                );

                            const fvPatchScalarField& rhop =
                                rho.boundaryField()[patchi];

                            forAll(rhop, fI)
                            {
                                rhoAdjustable += rhop[fI]*patchAreas[fI];
                            }
                        }
                    }
                    else
                    {
                        forAll(phip, i)
                        {
                            if (phip[i] < 0.0)
                            {
                                massIn -= phip[i];
                            }
                            else
                            {
                                fixedMassOut += phip[i];
                            }
                        }
                    }
                }
            }
        }

        reduce(areaAdjustable, sumOp<scalar>());
        if (areaAdjustable<SMALL) return false;

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = VSMALL + sum(mag(phi)).value();

        reduce(
            std::tie(massIn, fixedMassOut, adjustableMassOut),
            UniformParallelOp<sumOp<scalar>, 3>{}
        );

        if (compressible)
        {
            reduce(rhoAdjustable, sumOp<scalar>());
            rhoAdjustable /= areaAdjustable;
        }

        scalar massCorr = 0;
        scalar scaleCorr = 0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        if
        (
            magAdjustableMassOut > VSMALL
         && magAdjustableMassOut/totalFlux > SMALL
        )
        {
            scaleCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else
        {
            massCorr = massIn - fixedMassOut - adjustableMassOut;
        }

        forAll(bphi, patchi)
        {
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];

            if (!Up.coupled())
            {
                if
                (
                    isA<outflowVelocityFvPatchVectorField>(Up)
                )
                {
                    if (compressible)
                    {
                        massCorr /= (rhoAdjustable*areaAdjustable);
                    }
                    refCast<outflowVelocityFvPatchVectorField>
                    (
                        Up
                    ).scaleProfile(scaleCorr, massCorr);
                }
            }
        }
    }

    return outflowExist;
}

bool Foam::setOutFlowU
(
    const IOobject& rhoIO,
    surfaceScalarField& phi,
    volVectorField& U
)
{
    bool outflowExist = false;
    forAll(U.boundaryField(), pI)
    {
        const fvPatchVectorField& Up = U.boundaryField()[pI];
        if
        (
            isA<outflowVelocityFvPatchVectorField>(Up)
        )
        {
            outflowExist = true;
        }
    }

    if (outflowExist)
    {
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar areaAdjustable = 0.0;

        const surfaceScalarField::Boundary& bphi =
            phi.boundaryField();
        const fvMesh& mesh = phi.mesh();

        forAll(bphi, patchi)
        {
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled())
            {
                if
                (
                    Up.fixesValue()
                 && !isA<inletOutletFvPatchVectorField>(Up)
                 && !isA<interpolatedInletOutletFvPatchVectorField>(Up)
                 && !isA<windProfileDirectionVelocityFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }

                }
                else
                {
                    if
                    (
                        isA<outflowVelocityFvPatchVectorField>(Up)
                    )
                    {
                        areaAdjustable += gSum(mesh.boundary()[patchi].magSf());
                    }
                    else
                    {
                        forAll(phip, i)
                        {
                            if (phip[i] < 0.0)
                            {
                                massIn -= phip[i];
                            }
                            else
                            {
                                fixedMassOut += phip[i];
                            }
                        }
                    }
                }
            }
        }

        reduce(std::tie(massIn, fixedMassOut), UniformParallelOp<sumOp<scalar>, 2>{});

        scalar adjustableMassOut = massIn - fixedMassOut;

        forAll(bphi, patchi)
        {
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];

            if (!Up.coupled())
            {
                if
                (
                    isA<outflowVelocityFvPatchVectorField>(Up)
                )
                {
                    scalar patchArea = gSum(mesh.boundary()[patchi].magSf());
                    scalar thisPatchFlow =
                        patchArea/areaAdjustable*adjustableMassOut;

                    if (U.db().foundObject<volScalarField>(rhoIO.name()))
                    {
                        const volScalarField& rho =
                            U.db().lookupObject<volScalarField>
                            (
                                rhoIO.name()
                            );
                        scalar patchRho = gAverage(rho.boundaryField()[patchi]);
                        thisPatchFlow /= patchRho;
                    }
                    scalar averUn = thisPatchFlow/patchArea;

                    if (false)
                    {
                        Info<< "patch: " << mesh.boundary()[patchi].name()
                             << endl;
                        Info<< "Patch flow: " << thisPatchFlow << endl;
                        Info<< "Adjustable area: "  << areaAdjustable << endl;
                        Info<< "PatchArea: "  << patchArea << endl;
                        Info<< "Average velocity: "  << averUn << endl;
                    }

                    outflowVelocityFvPatchVectorField& Upatch =
                        refCast<outflowVelocityFvPatchVectorField>
                        (
                            Up
                        );
                    Upatch.forceAssign(mesh.boundary()[patchi].nf()*averUn);
                    phi.boundaryFieldRef()[patchi] =
                        Upatch&mesh.boundary()[patchi].Sf();
                }
            }
        }
    }

    return outflowExist;
}


// ************************************************************************* //
