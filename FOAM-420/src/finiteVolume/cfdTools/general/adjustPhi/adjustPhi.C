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
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/adjustPhi/adjustPhi.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchFields.H"
#include "fields/fvPatchFields/derived/windProfileDirectionVelocity/windProfileDirectionVelocityFvPatchVectorField.H"
#include "fields/fvPatchFields/derived/interpolatedInletOutlet/interpolatedInletOutletFvPatchFields.H"
#include "fields/fvPatchFields/basic/blended/blendedFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p
)
{
    if (p.needReference())
    {
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        surfaceScalarField::Boundary& bphi =
            phi.boundaryFieldRef();

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled())
            {
                if (isA<blendedFvPatchVectorField>(Up))
                {
                    List<bool> fixesValues(Up.fixesValues());

                    forAll(phip, i)
                    {
                        if (fixesValues[i])
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
                        else
                        {
                            if (phip[i] < 0.0)
                            {
                                massIn -= phip[i];
                            }
                            else
                            {
                                adjustableMassOut += phip[i];
                            }
                        }
                    }
                }
                else
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
                        forAll(phip, i)
                        {
                            if (phip[i] < 0.0)
                            {
                                massIn -= phip[i];
                            }
                            else
                            {
                                adjustableMassOut += phip[i];
                            }
                        }
                    }
                }
            }
        }

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = VSMALL + sum(mag(phi)).value();

        reduce(
            std::tie(massIn, fixedMassOut, adjustableMassOut),
            UniformParallelOp<sumOp<scalar>, 3>{}
        );

        scalar massCorr = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        if
        (
            magAdjustableMassOut > VSMALL
         && magAdjustableMassOut/totalFlux > SMALL
        )
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if (mag(fixedMassOut - massIn)/totalFlux > 1e-8)
        {
            WarningInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << endl;
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled())
            {
                if (isA<blendedFvPatchVectorField>(Up))
                {
                    List<bool> fixesValues(Up.fixesValues());
                    forAll(phip, i)
                    {
                        if (!fixesValues[i])
                        {
                            if (phip[i] > 0.0)
                            {
                                phip[i] *= massCorr;
                            }
                        }
                    }
                }
                else if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                }
            }
        }

        return mag(massIn)/totalFlux < SMALL
            && mag(fixedMassOut)/totalFlux < SMALL
            && mag(adjustableMassOut)/totalFlux < SMALL;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
