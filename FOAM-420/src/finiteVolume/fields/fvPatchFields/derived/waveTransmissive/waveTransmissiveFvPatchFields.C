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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/derived/waveTransmissive/waveTransmissiveFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/derived/waveTransmissive/waveTransmissiveFvPatchField.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(waveTransmissive);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::waveTransmissiveFvPatchField<Foam::scalar>::totalPressure()
{
    this->fieldInfMod_ = Field<scalar>(this->size(), 0.0);

    //hack to add total pressure-like behaveiour to this boundary
    if (totalPressure_)
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

        const fvsPatchField<scalar>& phip =
            patch().lookupPatchFieldInDb<surfaceScalarField, scalar>(db(), phiName_);

        const vectorField& Up
        (
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName_)
        );

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            this->fieldInfMod_ =
            (
                - 0.5*neg(phip)*magSqr(Up)
            );
        }
        else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            if (this->rhoName_ == "none")
            {
                const fvPatchField<scalar>& psip =
                    patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), psiName_);

                if (gamma_ > 1.0)
                {
                    scalar gM1ByG = (gamma_ - 1.0)/gamma_;

                    this->fieldInfMod_ =
                    (
                       fieldInf_
                       /pow
                        (
                            (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))
                            *magSqr(Up)),
                            1.0/gM1ByG
                        )
                    );
                }
                else
                {
                        this->fieldInfMod_ =
                        (
                            fieldInf_/(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up))
                        );
                }
            }
            else
            {
                const fvPatchField<scalar>& rho =
                    patch().lookupPatchFieldInDb<volScalarField, scalar>
                    (
                        db(),
                        this->rhoName_
                    );

                this->fieldInfMod_ =
                (
                    - 0.5*rho*neg(phip)*magSqr(Up)
                );
            }
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of " << phiName_ << " are not correct"
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
