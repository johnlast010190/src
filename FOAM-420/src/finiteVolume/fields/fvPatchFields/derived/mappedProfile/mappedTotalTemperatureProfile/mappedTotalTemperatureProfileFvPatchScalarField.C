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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedProfile/mappedTotalTemperatureProfile/mappedTotalTemperatureProfileFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "./fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchFields.H"


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

void Foam::mappedTotalTemperatureProfileFvPatchScalarField::updateValues()
{
    const fvPatchVectorField& Up =
        patch().lookupPatchFieldInDb<volVectorField, vector>
        (
            this->db(),
            UName_
        );

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    const fvPatchField<scalar>& psip =
        patch().lookupPatchFieldInDb<volScalarField, scalar>
        (
            this->db(),
            psiName_
        );

    scalar gM1ByG = (gamma_ - 1.0)/gamma_;

    forceAssign
    (
        this->mappedField()/
        (1.0 + 0.5*psip*gM1ByG*(neg(phip))*magSqr(Up))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedTotalTemperatureProfileFvPatchScalarField::
mappedTotalTemperatureProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedProfileFvPatchField<scalar>(p, iF),
    UName_("U"),
    phiName_("phi"),
    psiName_("thermo:psi"),
    gamma_(0.0)
{}


Foam::mappedTotalTemperatureProfileFvPatchScalarField::
mappedTotalTemperatureProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedProfileFvPatchField<scalar>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    gamma_(readScalar(dict.lookup("gamma")))
{}


Foam::mappedTotalTemperatureProfileFvPatchScalarField::
mappedTotalTemperatureProfileFvPatchScalarField
(
    const mappedTotalTemperatureProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedProfileFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_)
{}


Foam::mappedTotalTemperatureProfileFvPatchScalarField::
mappedTotalTemperatureProfileFvPatchScalarField
(
    const mappedTotalTemperatureProfileFvPatchScalarField& tppsf
)
:
    mappedProfileFvPatchField<scalar>(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_)
{}


Foam::mappedTotalTemperatureProfileFvPatchScalarField::
mappedTotalTemperatureProfileFvPatchScalarField
(
    const mappedTotalTemperatureProfileFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedProfileFvPatchField<scalar>(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedTotalTemperatureProfileFvPatchScalarField::write(Ostream& os)
const
{
    mappedProfileFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
    os.writeEntry("gamma", gamma_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedTotalTemperatureProfileFvPatchScalarField
    );
}

// ************************************************************************* //
