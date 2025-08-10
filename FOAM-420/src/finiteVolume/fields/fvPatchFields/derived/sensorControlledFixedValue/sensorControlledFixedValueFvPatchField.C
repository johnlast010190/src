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

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/sensorControlledFixedValue/sensorControlledFixedValueFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorControlledFixedValueFvPatchField<Type>::
sensorControlledFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    sensorName_("none"),
    bcValueController_(),
    isInlet_(true)
{}


template<class Type>
Foam::sensorControlledFixedValueFvPatchField<Type>::
sensorControlledFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    sensorName_(dict.lookup("sensor")),
    bcValueController_(Function1<Type>::New("Function", dict)),
    isInlet_(dict.lookupOrDefault<bool>("isInlet", true))
{
    //- initialise field
    if (dict.found("value"))
    {
        fvPatchField<Type>::forceAssign
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        Info<<"No value specified. Initialise from patchInternalField" << endl;
        fvPatchField<Type>::forceAssign(this->patchInternalField());
    }
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
Foam::sensorControlledFixedValueFvPatchField<Type>::
sensorControlledFixedValueFvPatchField
(
    const sensorControlledFixedValueFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    sensorName_(ptf.sensorName_),
    bcValueController_(ptf.bcValueController_, false),
    isInlet_(ptf.isInlet_)
{}


template<class Type>
Foam::sensorControlledFixedValueFvPatchField<Type>::
sensorControlledFixedValueFvPatchField
(
    const sensorControlledFixedValueFvPatchField& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    sensorName_(ptf.sensorName_),
    bcValueController_(ptf.bcValueController_, false),
    isInlet_(ptf.isInlet_)
{}


template<class Type>
Foam::sensorControlledFixedValueFvPatchField<Type>::
sensorControlledFixedValueFvPatchField
(
    const sensorControlledFixedValueFvPatchField& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    sensorName_(ptf.sensorName_),
    bcValueController_(),
    isInlet_(ptf.isInlet_)
{
    if (ptf.bcValueController_.valid())
    {
        bcValueController_.reset(ptf.bcValueController_->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::sensorControlledFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    if (!this->db().objectRegistry::template foundObject<sensor<scalar>>(sensorName_))
    {
        Info<<"creating sensor from dictionary" << endl;
        IOdictionary sensorDict
        (
            IOobject
            (
                "sensorDict",
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        autoPtr<sensor<scalar>> bcSensor = sensor<scalar>::New(mesh, sensorDict.subDict(sensorName_), sensorName_);
        regIOobject::store(bcSensor.ptr());
        Info<<"finished creating sensor" << endl;
    }

    sensor<scalar>& bcSensor_ = const_cast<sensor<scalar>&>
    (
        mesh.thisDb().lookupObject<sensor<scalar>>(sensorName_)
    );

    if (this->db().objectRegistry::template foundObject<volScalarField>(bcSensor_.fieldName()))
    {
        const volScalarField& sensorField(this->db().objectRegistry::template
            lookupObject<volScalarField>(bcSensor_.fieldName()));

        const scalar sensorVal = bcSensor_.value(sensorField);

        Type val = bcValueController_->value(sensorVal);

        // check sign of controller output field
        checkSign(val);

        fvPatchField<Type>::forceAssign(val);
        fixedValueFvPatchField<Type>::updateCoeffs();
    }
    else
    {
        Info<< "sensorControlledFixedValueFvPatchField::updateCoeffs()" << nl
            << "Cannot find volScalarField " << bcSensor_.fieldName()
            << " in mesh database!" << nl
            << "So this boundary condition never updates!" << endl;
    }
}


template<class Type>
void Foam::sensorControlledFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("sensor", sensorName_);
    if (bcValueController_.valid())
    {
        bcValueController_->writeData(os);
    }
    os.writeEntry("isInlet", isInlet_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
