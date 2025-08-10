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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "viscosityModels/approachingNu/approachingNu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(approachingNu, 0);
    addToRunTimeSelectionTable ( viscosityModel, approachingNu, dictionary );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
/* this is the linear blended algorithm... however I do prefer the alternative algorithm below
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::approachingNu::calcNu() const
{
    scalar current_time((U_.time()).value());
    Info<< "current_time = " << current_time << endl;
    Info<< "some_fixed_time_ = " << some_fixed_time_.value() << endl;

    if (current_time < some_fixed_time_.value())
    {
        scalar alpha = current_time / some_fixed_time_.value();
        Info<< "alpha = " << alpha << ", nu = " << (alpha * nu_desired_ + (1 - alpha) * nu_start_) << endl;
        return
        (
            strainRate() * dimensionedScalar("dummy",nu_desired_.dimensions()*dimTime,0.0) // dummy for adding strainRate dependend model
            + (alpha * nu_desired_ + (1 - alpha) * nu_start_)
        );
    }
    else
    {
        Info<< "nu = " << nu_desired_ << endl;
        return
        (
            strainRate() * dimensionedScalar("dummy",nu_desired_.dimensions()*dimTime,0.0) // dummy for adding strainRate dependend model
            + nu_desired_
        );
    }
}
*/
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::approachingNu::calcNu() const
{
    scalar current_time((U_.time()).value());

    if (current_time < some_fixed_time_.value())
    {
        scalar alpha
            = pow((nu_start_.value()*1000/nu_desired_.value()),
            (1/some_fixed_time_.value()));

        return
        (
            strainRate()
            * dimensionedScalar("dummy",nu_desired_.dimensions()*dimTime,0.0)
            + (nu_start_ / (pow(alpha, current_time))) + nu_desired_
        );
    }
    else
    {
        return
        (
            strainRate()
            * dimensionedScalar("dummy",nu_desired_.dimensions()*dimTime,0.0)
            + nu_desired_
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::approachingNu::approachingNu
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    approachingNuCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nu_desired_("nu_desired", dimViscosity, approachingNuCoeffs_),
    nu_start_("nu_start", dimViscosity, approachingNuCoeffs_),
    some_fixed_time_
    (
        "some_fixed_time", approachingNuCoeffs_.lookup("some_fixed_time")
    ),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::approachingNu::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    approachingNuCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    approachingNuCoeffs_.lookup("nu_desired") >> nu_desired_;
    approachingNuCoeffs_.lookup("nu_start") >> nu_start_;
    approachingNuCoeffs_.lookup("some_fixed_time") >> some_fixed_time_;

    return true;
}
