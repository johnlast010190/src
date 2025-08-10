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
    (c) 2016-2021 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "BokilBewtraMu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BokilBewtraMu, 0);
    addToRunTimeSelectionTable(materialModel, BokilBewtraMu, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BokilBewtraMu::BokilBewtraMu
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
    BokilBewtraMu::read();
    const word fName(this->phasePropertyName(XName_, phaseName));
    X_ = this->obr_.lookupObjectPtr<volScalarField>(fName);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::baseModels<Foam::scalar>* Foam::BokilBewtraMu::castScalarModel
(
    const word& modelName
)
{
    if (modelName == muModel::typeName)
    {
        return dynamic_cast<muModel*>(this);
    }
    return nullptr;
}


bool Foam::BokilBewtraMu::updateScalarField
(
    const word& fieldName,
    const volScalarField& volField
)
{
    if (fieldName == phasePropertyName(XName_, phaseName_))
    {
        X_ = &volField;
        return true;
    }
    return false;
}


Foam::tmp<Foam::scalarField> Foam::BokilBewtraMu::muPatch
(
    const label patchi
) const
{
    const scalarField XBc = X_->boundaryField()[patchi];
    if (X_ != nullptr)
    {
        const scalarField blending
        (
            min
            (
                max
                (
                    (Xc_ - XBc*fac_)/Xc_,
                    scalar(0.0)
                ),
                scalar(1.0)
            )
        );

        return min
        (
            mul_*blending
          + mu0_*pow(10.0, a_*XBc*fac_)*(1.0 - blending),
            mumax_
        );
    }

    return tmp<scalarField>(new scalarField(XBc.size(), mul_));
}


Foam::tmp<Foam::scalarField> Foam::BokilBewtraMu::muInternal() const
{
    const scalarField X = X_->primitiveField();
    if (X_ != nullptr)
    {
        const scalarField blending
        (
            min
            (
                max
                (
                    (Xc_ - X*fac_)/Xc_,
                    scalar(0.0)
                ),
                scalar(1.0)
            )
        );

        return min
        (
            mul_*blending
          + mu0_*pow(10.0, a_*X*fac_)*(1.0 - blending),
            mumax_
        );
    }

    return tmp<scalarField>(new scalarField(X.size(), mul_));
}


Foam::tmp<Foam::volScalarField> Foam::BokilBewtraMu::muGeometric() const
{
    tmp<volScalarField> tMu
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("mu", muModel::modelDims, Zero)
        )
    );
    volScalarField& mu = tMu.ref();

    mu.primitiveFieldRef() = muInternal();
    forAll(mu.boundaryField(), patchi)
    {
        mu.boundaryFieldRef()[patchi].forceAssign(muPatch(patchi));
    }

    return tMu;
}


Foam::scalar Foam::BokilBewtraMu::muCell(const label celli) const
{
    const scalar X = X_->operator[](celli);
    if (X_ != nullptr)
    {
        const scalar blending
        (
            min
            (
                max
                (
                    (Xc_ - X*fac_)/Xc_,
                    scalar(0.0)
                ),
                scalar(1.0)
            )
        );

        return min
        (
            mul_*blending
          + mu0_*pow(10.0, a_*X*fac_)*(1.0 - blending),
            mumax_
        );
    }

    return mul_;
}


bool Foam::BokilBewtraMu::read()
{
    XName_ = dict_->lookup<word>("XName");
    mu0_ = dict_->lookupOrDefault<scalar>("mu0", 0.00327);
    mul_ = dict_->lookupOrDefault<scalar>("mul", 1e-6);
    mumax_ = dict_->lookupOrDefault<scalar>("mumax", GREAT);
    Xc_ = dict_->lookupOrDefault<scalar>("Xc", 0.7);
    a_ = dict_->lookupOrDefault<scalar>("a", 0.132);
    fac_ = dict_->lookupOrDefault<scalar>("Xscaling", 1.0);

    return true;
}


// ************************************************************************* //
