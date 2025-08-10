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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "scalarModelFunction1.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/functions/Function1/Constant/Constant.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scalarModelFunction1, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        scalarModelFunction1,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarModelFunction1::scalarModelFunction1
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
    scalarModelFunction1::read();
}


Foam::autoPtr<Foam::scalarModelFunction1>
Foam::scalarModelFunction1::clone() const
{
    return autoPtr<scalarModelFunction1>
    (
        new scalarModelFunction1(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::scalarModelFunction1::castScalarModel
(
    const word& modelName
)
{
    return dynamic_cast<generalModel*>(this);
}


Foam::tmp<Foam::scalarField> Foam::scalarModelFunction1::generalInternal()
const
{
    if (useTime_)
    {
        return
            tmp<scalarField>
            (
                new scalarField
                (
                    mesh().nCells(),
                    general(mesh().time().value())
                )
            );
    }

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


Foam::tmp<Foam::scalarField> Foam::scalarModelFunction1::generalPatch
(
    const label patchi
) const
{
    if (useTime_)
    {
        return
            tmp<scalarField>
            (
                new scalarField
                (
                    mesh().boundary()[patchi].size(),
                    general(mesh().time().value())
                )
            );
    }

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


Foam::tmp<Foam::volScalarField> Foam::scalarModelFunction1::generalGeometric()
const
{
    if (useTime_)
    {
        return
            tmp<volScalarField>
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
                    dimensioned<scalar>
                    (
                        "general",
                        generalModel::dimensions(),
                        general(mesh().time().value())
                    )
                )
            );
    }

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


bool Foam::scalarModelFunction1::read()
{
    word dictName(dict_->dictName());
    fieldName_ = dict_->lookupOrDefault<word>("fieldName", "T");
    useTime_ = (fieldName_ == "time") ? true : false;
    if (!useTime_)
    {
        field_ = constructOrReturnRefFieldPtr<scalar>(fieldName_);
    }

    const word funcName(dictName.replace("ModelCoeffs", ""));
    value_.clear();
    value_ = Function1<scalar>::New(funcName, *dict_);

    generalModel::dimensions().reset
    (
        baseModels<scalar>::dimensions(funcName)
    );

    const word constType(Function1Types::Constant<scalar>::typeName);

    if
    (
        funcName == rhoModel::typeName
     && (value_->type() == constType || field_->isConst())
    )
    {
        isochoric_ = dict_->lookupOrDefault("isochoric", true);
    }
    else
    {
        isochoric_ = dict_->lookupOrDefault("isochoric", false);
    }

    if
    (
        funcName == rhoModel::typeName
     && (value_->type() == constType || fieldName_ != "p" || field_->isConst())
    )
    {
        incompressible_ = dict_->lookupOrDefault("incompressible", true);
    }
    else
    {
        incompressible_ = dict_->lookupOrDefault("incompressible", false);
    }

    // Having a physical property which has the same
    // value when measured in different directions.
    // (true always for scalar model)
    isotropic_ =  dict_->lookupOrDefault("isotropic", true);
    isConstant_ =
        dict_->lookupOrDefault
        (
            "isConstant",
            value_->type() == constType || field_->isConst()
        );

    return true;
}


// ************************************************************************* //
