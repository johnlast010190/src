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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "muSutherland.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(muSutherland, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        muSutherland,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dimensionedScalar Foam::muSutherland::readCoeff
(
    const word& coeffName,
    const dimensionSet& dims,
    const dictionary& dict
)
{
    return dimensionedScalar
    (
        coeffName,
        dims,
        readScalar(dict.subDict("transportCoeffs").lookup(coeffName))
    );
}


void Foam::muSutherland::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    const scalar rootT1 = sqrt(T1);
    const scalar mu1rootT2 = mu1*sqrt(T2);
    const scalar mu2rootT1 = mu2*rootT1;

    Ts_.value() = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_.value() = mu1*(1.0 + Ts_.value()/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::muSutherland::muSutherland
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    As_("As", dimDynamicViscosity/sqrt(dimTemperature), 0.0),
    Ts_("Ts", dimTemperature, 0.0)
{
    muSutherland::read();
}


Foam::autoPtr<Foam::muSutherland>
Foam::muSutherland::clone() const
{
    return autoPtr<muSutherland>
    (
        new muSutherland(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::muSutherland::castScalarModel
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


Foam::tmp<Foam::scalarField> Foam::muSutherland::muPatch
(
    const label patchi
) const
{
    const scalarField T(T_->boundaryField()[patchi]);
    return As_.value()*sqrt(T)/(1.0 + Ts_.value()/T);
}


Foam::tmp<Foam::scalarField> Foam::muSutherland::muInternal() const
{
    const scalarField T(T_->primitiveField());
    return As_.value()*sqrt(T)/(1.0 + Ts_.value()/T);
}


Foam::tmp<Foam::volScalarField> Foam::muSutherland::muGeometric() const
{
    volScalarField T(T_->operator()());
    return As_*sqrt(T)/(1.0 + Ts_/T);
}


Foam::scalar Foam::muSutherland::muCell(const label celli) const
{
    return
        As_.value()
       *::sqrt(T_->operator[](celli))
       /(1.0 + Ts_.value()/T_->operator[](celli));
}


bool Foam::muSutherland::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    As_.value() = dict_->lookup<scalar>("As");
    Ts_.value() = dict_->lookup<scalar>("Ts");

    return true;
}


// ************************************************************************* //
