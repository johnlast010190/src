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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "tabulated.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tabulated, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        tabulated,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulated::tabulated
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
    tabulated::read();
}


Foam::autoPtr<Foam::tabulated>
Foam::tabulated::clone() const
{
    return autoPtr<tabulated>
    (
        new tabulated(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::tabulated::castScalarModel
(
    const word& modelName
)
{
    return dynamic_cast<generalModel*>(this);
}


Foam::tmp<Foam::scalarField> Foam::tabulated::generalInternal() const
{
    tmp<scalarField> tInternalField(new scalarField(mesh().nCells()));
    scalarField& fieldCells = tInternalField.ref();

    // Initialize the internal fields
    tmp<scalarField> tfield1(field1_->primitiveField());
    tmp<scalarField> tfield2(field2_->primitiveField());
    scalarField& field1 = tfield1.ref();
    scalarField& field2 = tfield2.ref();

    forAll(fieldCells, celli)
    {
        fieldCells[celli] = general(field1[celli], field2[celli]);
    }
    return tInternalField;
}


Foam::tmp<Foam::scalarField> Foam::tabulated::generalPatch
(
    const label patchi
) const
{
    // Initialize the internal fields
    tmp<scalarField> tfield1(field1_->boundaryField()[patchi]);
    tmp<scalarField> tfield2(field2_->boundaryField()[patchi]);
    scalarField& field1 = tfield1.ref();
    scalarField& field2 = tfield2.ref();

    tmp<scalarField> tPatchField(new scalarField(field1.size()));
    scalarField& patchField = tPatchField.ref();

    forAll(patchField, facei)
    {
        patchField[facei] = general(field1[facei], field2[facei]);
    }
    return tPatchField;
}


Foam::tmp<Foam::volScalarField> Foam::tabulated::generalGeometric() const
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


bool Foam::tabulated::read()
{
    word dictName(dict_->dictName());
    fieldName1_ = dict_->lookupOrDefault<word>("fieldName1", "p");
    fieldName2_ = dict_->lookupOrDefault<word>("fieldName2", "T");
    field1_ = constructOrReturnRefFieldPtr<scalar>(fieldName1_);
    field2_ = constructOrReturnRefFieldPtr<scalar>(fieldName2_);

    const word funcName(dictName.replace("ModelCoeffs", ""));

    generalModel::dimensions().reset
    (
        baseModels<scalar>::dimensions(funcName)
    );

    table_.reset(new interpolation2DTable<scalar>(*dict_));

    incompressible_ = dict_->lookupOrDefault("incompressible", false);
    isochoric_ = dict_->lookupOrDefault("isochoric", false);
    isotropic_ = dict_->lookupOrDefault("isotropic", false);

    return true;
}


// ************************************************************************* //
