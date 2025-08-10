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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/basic/blended/blendedFvPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::blendedFvPatchField<Type>::blendedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    valueFraction_(p.size()),
    sensorName_("none"),
    bcValueController_(),
    updatedSwitch_(true),
    timeIndex_(0),
    boundaryOne_(),
    boundaryTwo_()
{}


template<class Type>
Foam::blendedFvPatchField<Type>::blendedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, true),
    valueFraction_("valueFraction", dict, p.size()),
    sensorName_(dict.lookup("sensor")),
    bcValueController_(),
    updatedSwitch_(false),
    timeIndex_(this->patch().boundaryMesh().mesh().time().timeIndex()),
    boundaryOne_(fvPatchField<Type>::New(p, iF, dict.subDict("boundaryOne")).ptr()),
    boundaryTwo_(fvPatchField<Type>::New(p, iF, dict.subDict("boundaryTwo")).ptr())
{
    if (dict.found("Function"))
    {
        bcValueController_ = Function1<scalar>::New("Function", dict);
    }

    // call evaluate causes problems e.g. in caseSetup if bc1 or bc2
    // are inletOutlet: phi cannot be found, then caseSetup crashes
    //evaluate();
}


template<class Type>
Foam::blendedFvPatchField<Type>::blendedFvPatchField
(
    const blendedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    valueFraction_(mapper(ptf.valueFraction_)),
    sensorName_(ptf.sensorName_),
    bcValueController_(ptf.bcValueController_, false),
    updatedSwitch_(ptf.updatedSwitch_),
    timeIndex_(ptf.timeIndex_),
    boundaryOne_(fvPatchField<Type>::New(ptf.boundaryOne(), p, iF, mapper).ptr()),
    boundaryTwo_(fvPatchField<Type>::New(ptf.boundaryTwo(), p, iF, mapper).ptr())
{
    if (!mapper.suppressWarning() && notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::blendedFvPatchField<Type>::blendedFvPatchField
(
    const blendedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    valueFraction_(ptf.valueFraction_),
    sensorName_(ptf.sensorName_),
    bcValueController_(ptf.bcValueController_, false),
    updatedSwitch_(ptf.updatedSwitch_),
    timeIndex_(ptf.timeIndex_),
    boundaryOne_(ptf.boundaryOne().clone().ptr()),
    boundaryTwo_(ptf.boundaryTwo().clone().ptr())
{}


template<class Type>
Foam::blendedFvPatchField<Type>::blendedFvPatchField
(
    const blendedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    valueFraction_(ptf.valueFraction_),
    sensorName_(ptf.sensorName_),
    bcValueController_(),
    updatedSwitch_(ptf.updatedSwitch_),
    timeIndex_(ptf.timeIndex_),
    boundaryOne_(ptf.boundaryOne().clone().ptr()),
    boundaryTwo_(ptf.boundaryTwo().clone().ptr())
{
    if (ptf.bcValueController_.valid())
    {
        bcValueController_.reset(ptf.bcValueController_->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::blendedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (debug)
    {
        Info<<"call to Foam::blendedFvPatchField<Type>::autoMap(m)" << endl;
    }

    boundaryOne_().autoMap(m);
    boundaryTwo_().autoMap(m);

    fvPatchField<Type>::autoMap(m);
    m(valueFraction_, valueFraction_);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    if (debug)
    {
        Info<<"call to Foam::blendedFvPatchField<Type>::rmap(ptf, addr)" << endl;
    }

    boundaryOne_().rmap(refCast<const blendedFvPatchField<Type>>(ptf).boundaryOne_(), addr);
    boundaryTwo_().rmap(refCast<const blendedFvPatchField<Type>>(ptf).boundaryTwo_(), addr);

    fvPatchField<Type>::rmap(ptf, addr);

    const blendedFvPatchField<Type>& mptf =
        refCast<const blendedFvPatchField<Type>>(ptf);

    valueFraction_.rmap(mptf.valueFraction_, addr);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fvPatchField<Type>::autoMapGIB(mapper);

    boundaryOne_().autoMapGIB(mapper);
    boundaryTwo_().autoMapGIB(mapper);

    mapper.map(valueFraction_, scalar(0));
}


template<class Type>
void Foam::blendedFvPatchField<Type>::calcSwitch()
{
    // do not update the switch at constructor call
    // reason: when this boundary condition is constructed for the field
    // which is used as switching criterion, evaluation does not work
    if (!constructorCall_)
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();

        // update only once per time step
        if (timeIndex_ != mesh.time().timeIndex())
        {
            updatedSwitch_ = false;
            timeIndex_ = mesh.time().timeIndex();
        }

        // instantiate or lookup sensor:
        //if (!mesh.thisDb().foundObject<sensor<scalar>>(sensorName_))
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

        if (!updatedSwitch_)
        {
            sensor<scalar>& bcSensor_ =
                mesh.thisDb().lookupObjectRef<sensor<scalar>>(sensorName_);

            if (this->db().objectRegistry::template foundObject<volScalarField>(bcSensor_.fieldName()))
            {
                if (debug)
                {
                    Info<<"calculating switch for field : " << this->internalField().name() << endl;
                }

                const scalarField sensorField(bcSensor_.valueField());

                if (sensorField.size() != this->patch().size())
                {
                    FatalErrorInFunction
                        << "Sensor '" << sensorName_
                        << "' returned inconsistent size for patch '" << this->patch().name()
                        << "'" << exit(FatalError);
                }

                forAll(sensorField, i)
                {
                    valueFraction_[i] = bcValueController_->value(sensorField[i]);
                }

                updatedSwitch_ = true;
                timeIndex_ = this->patch().boundaryMesh().mesh().time().timeIndex();
            }
        }
    }
}


template<class Type>
void Foam::blendedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (debug)
    {
        Info<<"call Foam::blendedFvPatchField<Type>::evaluate(commsTypes)" << endl;
    }

    boundaryOne_().evaluate();
    boundaryTwo_().evaluate();

    calcSwitch();

    if (!this->updated())
    {
        updateCoeffs();
    }

    Field<Type>::operator=
    (
        boundaryOne_() * valueFraction_
      + boundaryTwo_() * (1. - valueFraction_)
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::blendedFvPatchField<Type>::updateCoeffs()
{
    if (debug)
    {
        Info<<"call Foam::blendedFvPatchField<Type>::updateCoeffs()" << endl;
    }

    boundaryOne_().updateCoeffs();
    boundaryTwo_().updateCoeffs();

    calcSwitch();

    Field<Type>::operator=
    (
        boundaryOne_() * valueFraction_
      + boundaryTwo_() * (1. - valueFraction_)
    );

    constructorCall_ = false;
    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::blendedFvPatchField<Type>::snGrad() const
{
    return boundaryOne_().snGrad() * valueFraction_
         + boundaryTwo_().snGrad() * (1. - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::blendedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& tsf
) const
{
    return boundaryOne_().valueInternalCoeffs(tsf) * valueFraction_
         + boundaryTwo_().valueInternalCoeffs(tsf) * (1. - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::blendedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& tsf
) const
{
    return boundaryOne_().valueBoundaryCoeffs(tsf) * valueFraction_
         + boundaryTwo_().valueBoundaryCoeffs(tsf) * (1. - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::blendedFvPatchField<Type>::gradientInternalCoeffs() const
{
    return boundaryOne_().gradientInternalCoeffs() * valueFraction_
         + boundaryTwo_().gradientInternalCoeffs() * (1. - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::blendedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return boundaryOne_().gradientBoundaryCoeffs() * valueFraction_
         + boundaryTwo_().gradientBoundaryCoeffs() * (1. - valueFraction_);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("sensor", sensorName_);
    if (bcValueController_.valid())
    {
        bcValueController_->writeData(os);
    }

    // Write bc 1
    os.beginBlock("boundaryOne");
    boundaryOne_().write(os);
    os.endBlock();

    // Write bc 2
    os.beginBlock("boundaryTwo");
    boundaryTwo_().write(os);
    os.endBlock();

    // Write members
    valueFraction_.writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::blendedFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    boundaryOne_().operator=(ul);
    boundaryTwo_().operator=(ul);
    fvPatchField<Type>::operator=(ul);
}

template<class Type>
void Foam::blendedFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    boundaryOne_().operator=(ptf);
    boundaryTwo_().operator=(ptf);
    fvPatchField<Type>::operator=(ptf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    boundaryOne_().operator+=(ptf*valueFraction_);
    boundaryTwo_().operator+=(ptf*(1.-valueFraction_));
    fvPatchField<Type>::operator+=(ptf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    boundaryOne_().operator-=(ptf*valueFraction_);
    boundaryTwo_().operator-=(ptf*(1.-valueFraction_));
    fvPatchField<Type>::operator-=(ptf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    boundaryOne_().operator*=(ptf);
    boundaryTwo_().operator*=(ptf);
    fvPatchField<Type>::operator*=(ptf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    boundaryOne_().operator/=(ptf);
    boundaryTwo_().operator/=(ptf);
    fvPatchField<Type>::operator/=(ptf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    boundaryOne_().operator+=(tf*valueFraction_);
    boundaryTwo_().operator+=(tf*(1.-valueFraction_));
    fvPatchField<Type>::operator+=(tf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    boundaryOne_().operator-=(tf*valueFraction_);
    boundaryTwo_().operator-=(tf*(1.-valueFraction_));
    fvPatchField<Type>::operator-=(tf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    boundaryOne_().operator*=(tf);
    boundaryTwo_().operator*=(tf);
    fvPatchField<Type>::operator*=(tf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    boundaryOne_().operator/=(tf);
    boundaryTwo_().operator/=(tf);
    fvPatchField<Type>::operator/=(tf);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator=
(
    const Type& t
)
{
    boundaryOne_().operator=(t);
    boundaryTwo_().operator=(t);
    fvPatchField<Type>::operator=(t);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator+=
(
    const Type& t
)
{
    boundaryOne_().operator+=(t*valueFraction_);
    boundaryTwo_().operator+=(t*(1.-valueFraction_));
    fvPatchField<Type>::operator+=(t);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator-=
(
    const Type& t
)
{
    boundaryOne_().operator-=(t*valueFraction_);
    boundaryTwo_().operator-=(t*(1.-valueFraction_));
    fvPatchField<Type>::operator-=(t);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator*=
(
    const scalar s
)
{
    boundaryOne_().operator*=(s);
    boundaryTwo_().operator*=(s);
    fvPatchField<Type>::operator*=(s);
}


template<class Type>
void Foam::blendedFvPatchField<Type>::operator/=
(
    const scalar s
)
{
    boundaryOne_().operator/=(s);
    boundaryTwo_().operator/=(s);
    fvPatchField<Type>::operator/=(s);
}

// ************************************************************************* //
