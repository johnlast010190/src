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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/adjustSplitBCs/adjustSplitBCs.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchFields.H"
#include "fields/fvPatchFields/derived/windProfileDirectionVelocity/windProfileDirectionVelocityFvPatchVectorField.H"
#include "fields/fvPatchFields/derived/targetFlowRateOutletPressure/targetFlowRateOutletPressureFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/interpolatedInletOutlet/interpolatedInletOutletFvPatchFields.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::adjustSplitBCs
(
    const surfaceScalarField& phi,
    const volVectorField& U,
    GeometricField<Type, fvPatchField, volMesh>& Up
)
{
    bool massFlowSplitExist = false;
    forAll(Up.boundaryField(), pI)
    {
        const fvPatchField<Type>& Upp = Up.boundaryField()[pI];
        if (isA<flowRateBase>(Upp))
        {
            if (refCast<const flowRateBase>(Upp).isSplitBC())
            {
                massFlowSplitExist = true;
            }
        }
    }

    if (massFlowSplitExist)
    {
        const surfaceScalarField::Boundary& bphi = phi.boundaryField();

        // Fixed inflow mass per boundary patch
        List<scalar> massIn(bphi.size(), 0.0);

        label lastUncoupledPatch = 0;
        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled())
            {
                lastUncoupledPatch = patchi;
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
                            massIn[patchi] -= phip[i];
                        }
                    }

                }
            }
        }

        // bphi.size() will not be the same on all processors due to processor
        // patches, so we shrink it
        massIn.resize(lastUncoupledPatch+1, 0.0);
        listReduce(massIn, sumOp<scalar>());

        forAll(bphi, patchi)
        {
            fvPatchField<Type>& Upp = Up.boundaryFieldRef()[patchi];

            if (!Upp.coupled())
            {
                if (isA<flowRateBase>(Upp))
                {
                    flowRateBase& tUpp = refCast<flowRateBase>(Upp);
                    if (tUpp.isSplitBC())
                    {
                        tUpp.updateInletFlowRate(massIn);
                    }
                }
            }
        }
    }

    return massFlowSplitExist;
}


// ************************************************************************* //
