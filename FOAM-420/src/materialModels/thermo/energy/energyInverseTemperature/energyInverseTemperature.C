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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "energyInverseTemperature.H"
#include "thermo/energy/matAbsoluteEnthalpy/matAbsoluteEnthalpy.H"
#include "thermo/energy/matAbsoluteInternalEnergy/matAbsoluteInternalEnergy.H"
#include "thermo/energy/matSensibleEnthalpy/matSensibleEnthalpy.H"
#include "thermo/energy/matSensibleInternalEnergy/matSensibleInternalEnergy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * //

const Foam::scalar Foam::energyInverseTemperature::tol_ = 1.0e-4;

const int Foam::energyInverseTemperature::maxIter_ = 100;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(energyInverseTemperature, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        energyInverseTemperature,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::energyInverseTemperature::lookupEnergyFieldName()
{
    const word eName = materialsDict().lookup<word>("energy");
    if (eName == matSensibleInternalEnergy::typeName)
    {
        return matSensibleInternalEnergy::name();
    }
    else if (eName == matSensibleEnthalpy::typeName)
    {
        return matSensibleEnthalpy::name();
    }
    else if (eName == matAbsoluteInternalEnergy::typeName)
    {
        return matAbsoluteInternalEnergy::name();
    }
    else if (eName == matAbsoluteEnthalpy::typeName)
    {
        return matAbsoluteEnthalpy::name();
    }
    else
    {
        FatalErrorInFunction
            << "Couldn't guess energy type for "
            << "the temperature inverse function. " << nl
            << exit(FatalError);
    }

    return word::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyInverseTemperature::energyInverseTemperature
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
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(5);
    energyInverseTemperature::read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::energyInverseTemperature::~energyInverseTemperature()
{}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

void Foam::energyInverseTemperature::updateTable(const word& modelName)
{
    const word mixtureName =
        materialTables_.tableName(phaseName_, specieName_);
    if (!materialTables_.sTable().found(mixtureName))
    {
        FatalErrorInFunction
            << "Temperature model requires: \"" << mixtureName
            << "\" model."
            << abort(FatalError);
    }

    const matScalarTable& models = materialTables_.sTable()[mixtureName];

    if
    (
        modelName == TModel::typeName
     && models.found(HEModel::typeName)
     && models.found(CpvModel::typeName)
    )
    {
        sMod_.set(HE, models[HEModel::typeName]);
        sMod_.set(Cpv, models[CpvModel::typeName]);
        dep_[0].model = models[TModel::typeName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, models[HEModel::typeName]);
        dep_[0].dependencies.set(1, models[CpvModel::typeName]);

        // Limit model is optional
        if (models.found(limitModel::typeName))
        {
            sMod_.set(limit, models[limitModel::typeName]);
        }
    }
    else if
    (
        modelName == THaModel::typeName
     && models.found(HaModel::typeName)
     && models.found(CpModel::typeName)
    )
    {
        sMod_.set(Ha, models[HaModel::typeName]);
        sMod_.set(Cp, models[CpModel::typeName]);
        dep_[1].model = models[THaModel::typeName];
        dep_[1].dependencies.setSize(2);
        dep_[1].dependencies.set(0, models[HaModel::typeName]);
        dep_[1].dependencies.set(1, models[CpModel::typeName]);

        // Limit model is optional
        if (models.found(limitModel::typeName))
        {
            sMod_.set(limit, models[limitModel::typeName]);
        }
    }
    else if
    (
        modelName == THsModel::typeName
     && models.found(HsModel::typeName)
     && models.found(CpModel::typeName)
    )
    {
        sMod_.set(Hs, models[HsModel::typeName]);
        sMod_.set(Cp, models[CpModel::typeName]);
        dep_[2].model = models[THsModel::typeName];
        dep_[2].dependencies.setSize(2);
        dep_[2].dependencies.set(0, models[HsModel::typeName]);
        dep_[2].dependencies.set(1, models[CpModel::typeName]);

        // Limit model is optional
        if (models.found(limitModel::typeName))
        {
            sMod_.set(limit, models[limitModel::typeName]);
        }
    }
    else if
    (
        modelName == TEaModel::typeName
     && models.found(EaModel::typeName)
     && models.found(CvModel::typeName)
    )
    {
        sMod_.set(Ea, models[EaModel::typeName]);
        sMod_.set(Cv, models[CvModel::typeName]);
        dep_[3].model = models[TEaModel::typeName];
        dep_[3].dependencies.setSize(2);
        dep_[3].dependencies.set(0, models[EaModel::typeName]);
        dep_[3].dependencies.set(1, models[CvModel::typeName]);

        // Limit model is optional
        if (models.found(limitModel::typeName))
        {
            sMod_.set(limit, models[limitModel::typeName]);
        }
    }
    else if
    (
        modelName == TEsModel::typeName
     && models.found(EsModel::typeName)
     && models.found(CvModel::typeName)
    )
    {
        sMod_.set(Es, models[EsModel::typeName]);
        sMod_.set(Cv, models[CvModel::typeName]);
        dep_[4].model = models[TEsModel::typeName];
        dep_[4].dependencies.setSize(2);
        dep_[4].dependencies.set(0, models[EsModel::typeName]);
        dep_[4].dependencies.set(1, models[CvModel::typeName]);

        // Limit model is optional
        if (models.found(limitModel::typeName))
        {
            sMod_.set(limit, models[limitModel::typeName]);
        }
    }
    else
    {
        FatalErrorInFunction
            << modelName << " model requires: " << mixtureName
            << " " << HEModel::typeName << " and " << CpvModel::typeName << nl
            << " available models are: " << models.toc() << nl
            << " in table: \"" << mixtureName << "\"" << nl
            << abort(FatalError);
    }
}


Foam::baseModels<Foam::scalar>* Foam::energyInverseTemperature::castScalarModel
(
    const word& modelName
)
{
    if (modelName == TModel::typeName)
    {
        return dynamic_cast<TModel*>(this);
    }
    else if (modelName == THsModel::typeName)
    {
        return dynamic_cast<THsModel*>(this);
    }
    else if (modelName == THaModel::typeName)
    {
        return dynamic_cast<THaModel*>(this);
    }
    else if (modelName == TEaModel::typeName)
    {
        return dynamic_cast<TEaModel*>(this);
    }
    else if (modelName == TEsModel::typeName)
    {
        return dynamic_cast<TEsModel*>(this);
    }
    return nullptr;
}


bool Foam::energyInverseTemperature::updateScalarField
(
    const word& fieldName,
    const volScalarField& volField
)
{
    if (fieldName == phasePropertyName("T"))
    {
        TnonConst_ = const_cast<volScalarField*>(&volField);
        return true;
    }
    else if (fieldName == phasePropertyName(lookupEnergyFieldName()))
    {
        he_ = &volField;
        return true;
    }
    return false;
}


Foam::scalar Foam::energyInverseTemperature::Cell
(
    const label celli,
    const label HEFun,
    const label CpvFun
) const
{
    scalar& Test = TnonConst_->operator[](celli);
    if ((Test + Tref_) < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: "
            << (Test + Tref_) << " [K]."
            << abort(FatalError);
    }
    else if (T_->isConst())
    {
        return Tref_ + T_->constant().value();
    }

    //- Restore temperature after use
    const scalar Trestore = Test;

    scalar Tnew = Test;
    scalar Ttol = (Test + Tref_)*tol_;
    int iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (sMod_[HEFun][celli] - he_->operator[](celli))
           /sMod_[CpvFun][celli];

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    if (sMod_(limit) != nullptr)
    {
        Tnew = sMod_[limit][celli];
    }

    //- Restore temperature in registry
    Test = Trestore;

    return Tnew;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::Patch
(
    const label patchi,
    const label HEFun,
    const label CpvFun
) const
{
    scalarField& Test = TnonConst_->boundaryFieldRef()[patchi];
    if (TnonConst_->boundaryFieldRef()[patchi].empty())
    {
        return tmp<scalarField>(new scalarField());
    }
    const scalarField& he = he_->boundaryField()[patchi];
    if ((min(Test) + Tref_) < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: "
            << (min(Test) + Tref_) << " [K]."
            << abort(FatalError);
    }
    else if (T_->isConst())
    {
        return
            tmp<scalarField>
            (
                new scalarField(Test.size(), Tref_ + T_->constant().value())
            );
    }

    // Store old value of temperature
    scalarField Told(Test);

    tmp<scalarField> tTnew(new scalarField(Test));
    scalarField& Tnew = tTnew.ref();
    const scalar Ttol((min(Test) + Tref_)*tol_);
    int iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (sMod_[HEFun].boundaryField()[patchi] - he)
           /sMod_[CpvFun].boundaryField()[patchi];

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }
    } while (max(mag(Tnew - Test)) > Ttol);

    if (sMod_(limit) != nullptr)
    {
        Tnew = sMod_[limit].boundaryField()[patchi];
    }

    // Restore temperature in registry
    Test = Told;
    return tTnew;
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::Internal
(
    const label HEFun,
    const label CpvFun
) const
{
    scalarField& Test = TnonConst_->primitiveFieldRef();
    if (!Test.size())
    {
        return tmp<scalarField>(new scalarField());
    }
    const scalarField& he = he_->primitiveField();
    if ((min(Test) + Tref_) < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: "
            << (min(Test) + Tref_) << " [K]."
            << abort(FatalError);
    }
    else if (T_->isConst())
    {
        return
            tmp<scalarField>
            (
                new scalarField(Test.size(), Tref_ + T_->constant().value())
            );
    }

    // Store old value of temperature
    scalarField Told(Test);

    tmp<scalarField> tTnew(new scalarField(Test));
    scalarField& Tnew = tTnew.ref();
    const scalar Ttol((min(Test) + Tref_)*tol_);

    int iter = 0;
    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (sMod_[HEFun].primitiveField() - he)
           /sMod_[CpvFun].primitiveField();
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }
    } while (max(mag(Tnew - Test)) > Ttol);

    if (sMod_(limit) != nullptr)
    {
        Tnew = sMod_[limit].primitiveField();
    }

    // Restore temperature in registry
    Test = Told;
    return tTnew;
}


Foam::scalar Foam::energyInverseTemperature::TCell(const label celli) const
{
    return Cell(celli, HE, Cpv);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TPatch
(
    const label patchi
) const
{
    return Patch(patchi, HE, Cpv);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TInternal() const
{
    return Internal(HE, Cpv);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::TGeometric() const
{
    tmp<volScalarField> tTest(new volScalarField("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(HE, Cpv);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, HE, Cpv));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::THaCell(const label celli) const
{
    return Cell(celli, Ha, Cp);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::THaPatch
(
    const label patchi
) const
{
    return Patch(patchi, Ha, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::THaInternal() const
{
    return Internal(Ha, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::THaGeometric() const
{
    tmp<volScalarField> tTest(new volScalarField("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(Ha, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, Ha, Cp));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::THsCell(const label celli) const
{
    return Cell(celli, Hs, Cp);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::THsPatch
(
    const label patchi
) const
{
    return Patch(patchi, Hs, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::THsInternal() const
{
    return Internal(Hs, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::THsGeometric() const
{
    tmp<volScalarField> tTest(new volScalarField("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(Hs, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, Hs, Cp));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::TEaCell(const label celli) const
{
    return Cell(celli, Ea, Cp);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TEaPatch
(
    const label patchi
) const
{
    return Patch(patchi, Ea, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::TEaInternal() const
{
    return Internal(Ea, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::TEaGeometric() const
{
    tmp<volScalarField> tTest(new volScalarField("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(Ea, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, Ea, Cp));
    }
    return tTest;
}


Foam::scalar Foam::energyInverseTemperature::TEsCell(const label celli) const
{
    return Cell(celli, Es, Cp);
}


Foam::tmp<Foam::scalarField> Foam::energyInverseTemperature::TEsPatch
(
    const label patchi
) const
{
    return Patch(patchi, Es, Cp);
}


Foam::tmp<Foam::scalarField>
Foam::energyInverseTemperature::TEsInternal() const
{
    return Internal(Es, Cp);
}


Foam::tmp<Foam::volScalarField>
Foam::energyInverseTemperature::TEsGeometric() const
{
    tmp<volScalarField> tTest(new volScalarField("Test", *TnonConst_));
    volScalarField& Test = tTest.ref();
    Test.primitiveFieldRef() = Internal(Es, Cp);
    forAll(Test.boundaryField(), patchi)
    {
        Test.boundaryFieldRef()[patchi].forceAssign(Patch(patchi, Es, Cp));
    }
    return tTest;
}


bool Foam::energyInverseTemperature::read()
{
    TnonConst_ = lookupPtr<scalar>("T");
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    Tref_ = T_->offset().value();
    he_ = lookupPtr<scalar>(lookupEnergyFieldName());

    return true;
}


// ************************************************************************* //
