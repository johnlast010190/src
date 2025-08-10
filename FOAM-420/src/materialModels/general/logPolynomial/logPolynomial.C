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
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "logPolynomial.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(logPolynomial, 0);
    addToRunTimeSelectionTable(materialModel, logPolynomial, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::logPolynomial::logPolynomial
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    logPolynomial::read();
}


Foam::autoPtr<Foam::logPolynomial>
Foam::logPolynomial::clone() const
{
    return autoPtr<logPolynomial>
    (
        new logPolynomial(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::logPolynomial::castScalarModel
(
    const word& modelName
)
{
    return dynamic_cast<generalModel*>(this);
}


Foam::tmp<Foam::scalarField> Foam::logPolynomial::generalInternal() const
{
    tmp<scalarField> tInternalField(new scalarField(mesh().nCells()));
    scalarField& fieldCells = tInternalField.ref();
    tmp<scalarField> tfield(field_->primitiveField());
    scalarField& field = tfield.ref();
    forAll(fieldCells, celli)
    {
        fieldCells[celli] = general(field[celli]);
    }
    return tInternalField;
}


Foam::tmp<Foam::scalarField> Foam::logPolynomial::generalPatch
(
    const label patchi
) const
{
    tmp<scalarField> tfield(field_->boundaryField()[patchi]);
    scalarField& field = tfield.ref();
    tmp<scalarField> tPatchField(new scalarField(field.size()));
    scalarField& patchField = tPatchField.ref();
    forAll(patchField, facei)
    {
        patchField[facei] = general(field[facei]);
    }
    return tPatchField;
}


Foam::tmp<Foam::volScalarField> Foam::logPolynomial::generalGeometric() const
{
    tmp<volScalarField> tField
    (
        new volScalarField
        (
            IOobject
            (
                "general",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensioned<scalar>("general", generalModel::dimensions(), Zero)
        )
    );
    volScalarField& field = tField.ref();
    field.primitiveFieldRef() = generalInternal();
    forAll(field.boundaryField(), patchi)
    {
        field.boundaryFieldRef()[patchi].forceAssign(generalPatch(patchi));
    }
    return tField;
}


bool Foam::logPolynomial::read()
{
    word dictName(dict_->dictName());
    fieldName_ = dict_->lookupOrDefault<word>("fieldName", "T");
    field_ = constructOrReturnRefFieldPtr<scalar>(fieldName_);

    const word funcName(dictName.replace("ModelCoeffs", ""));
    coeffs_.reset(dict_->lookup<scalarField>(funcName + "Coeffs"));

    generalModel::dimensions().reset
    (
        baseModels<scalar>::dimensions(funcName)
    );

    incompressible_ = dict_->lookupOrDefault("incompressible", false);
    isochoric_ = dict_->lookupOrDefault("isochoric", false);
    isotropic_ = dict_->lookupOrDefault("isotropic", false);

    return true;
}


// ************************************************************************* //
