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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/waveTransmissive/waveTransmissiveFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/CrankNicolsonDdtScheme/CrankNicolsonDdtScheme.H"
#include "finiteVolume/ddtSchemes/backwardDdtScheme/backwardDdtScheme.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::waveTransmissiveFvPatchField<Type>::totalPressure()
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("thermo:psi"),
    gamma_(0.0),
    totalPressure_(false),
    UName_("U")
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    totalPressure_(ptf.totalPressure_),
    UName_(ptf.UName_)
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    totalPressure_(dict.lookupOrDefault<bool>("totalPressure", false)),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_),
    totalPressure_(ptpsf.totalPressure_),
    UName_(ptpsf.UName_)
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_),
    totalPressure_(ptpsf.totalPressure_),
    UName_(ptpsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::waveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip =
        this->patch().template
            lookupPatchFieldInDb<volScalarField, scalar>(this->db(), psiName_);

    const surfaceScalarField& phi =
        this->db().template lookupObject<surfaceScalarField>(this->phiName_);

    fvsPatchField<scalar> phip =
        this->patch().template
            lookupPatchFieldInDb<surfaceScalarField, scalar>(this->db(), this->phiName_);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop =
            this->patch().template
                lookupPatchFieldInDb<volScalarField, scalar>(this->db(), this->rhoName_);

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf() + sqrt(gamma_/psip);
}

template<class Type>
void Foam::waveTransmissiveFvPatchField<Type>::updateCoeffs()
{
    if
    (
        this->db().objectRegistry::template foundObject<surfaceScalarField>
        (advectiveFvPatchField<Type>::phiName_)
     && this->db().objectRegistry::template foundObject<volScalarField>
        (advectiveFvPatchField<Type>::rhoName_)
     && this->db().objectRegistry::template foundObject<volScalarField>
        (psiName_)
    )
    {
        if (totalPressure_) totalPressure();

        advectiveFvPatchField<Type>::updateCoeffs();
    }
    else
    {
        this->valueFraction() = 1.0;
        mixedFvPatchField<Type>::updateCoeffs();
    }
}

template<class Type>
void Foam::waveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    this->template
        writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    this->template
        writeEntryIfDifferent<word>(os, "rho", "rho", this->rhoName_);
    this->template
        writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
    this->template
        writeEntryIfDifferent<word>(os, "U", "U", UName_);
    this->template
        writeEntryIfDifferent<bool>(os, "totalPressure", false, totalPressure_);

    os.writeEntry("gamma", gamma_);

    if (this->lInf_ > SMALL)
    {
        os.writeEntry("fieldInf", this->fieldInf_);
        os.writeEntry("lInf",this->lInf_);
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
