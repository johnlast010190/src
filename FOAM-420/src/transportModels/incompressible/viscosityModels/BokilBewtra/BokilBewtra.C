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
    (c) 2016 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "viscosityModels/BokilBewtra/BokilBewtra.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(BokilBewtra, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        BokilBewtra,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::BokilBewtra::calcNu() const
{
    if (U_.db().foundObject<volScalarField>(XName_))
    {
        // get solid tracer field
        const volScalarField& X = U_.db().lookupObject<volScalarField>(XName_);

        // make sure exponent doesnt have dimensions but value is in mg/L
        dimensionedScalar fac("fac", Xc_.dimensions()/X.dimensions(), fac_);

        volScalarField blending
        (
            min
            (
                max((Xc_-X*fac)/Xc_, scalar(0.0)), scalar(1.0)
            )
        );

        return min
        (
            nul_ * blending
          + nu0_ * pow(10.0, a_*X*fac) * (1.-blending),
            numax_
        );
    }

    Info<< "Field " << XName_ << " not found. Return liquid viscosity." << endl;
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nul",
                U_.time().timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            nul_
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::BokilBewtra::BokilBewtra
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    BokilBewtraCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    XName_(BokilBewtraCoeffs_.lookup("XName")),
    nu0_
    (
        dimensionedScalar::lookupOrDefault
        (
            "nu0",
            BokilBewtraCoeffs_,
            dimViscosity,
            3.27e-6
        )
    ),
    nul_
    (
        dimensionedScalar::lookupOrDefault
        (
            "nul",
            BokilBewtraCoeffs_,
            dimViscosity,
            1e-6
        )
    ),
    numax_
    (
        dimensionedScalar::lookupOrDefault
        (
            "numax",
            BokilBewtraCoeffs_,
            dimViscosity,
            GREAT
        )
    ),
    Xc_
    (
        dimensionedScalar::lookupOrDefault
        (
            "Xc",
            BokilBewtraCoeffs_,
            dimDensity,
            0.7
        )
    ),
    a_
    (
        dimensionedScalar::lookupOrDefault
        (
            "a",
            BokilBewtraCoeffs_,
            dimless/Xc_.dimensions(),
            0.132
        )
    ),
    fac_(BokilBewtraCoeffs_.lookupOrDefault<scalar>("Xscaling", 1.0)),
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
        //calcNu() // X not yet available at construction!
        U_.mesh(),
        nul_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::BokilBewtra::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    BokilBewtraCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    XName_ = word(BokilBewtraCoeffs_.lookup("XName"));

    nu0_ = dimensionedScalar::lookupOrDefault
    (
        "nu0",
        BokilBewtraCoeffs_,
        dimViscosity,
        0.00327
    );

    nul_ = dimensionedScalar::lookupOrDefault
    (
        "nul",
        BokilBewtraCoeffs_,
        dimViscosity,
        1e-6
    );

    numax_ = dimensionedScalar::lookupOrDefault
    (
        "numax",
        BokilBewtraCoeffs_,
        dimViscosity,
        GREAT
    );

    Xc_ = dimensionedScalar::lookupOrDefault
    (
        "Xc",
        BokilBewtraCoeffs_,
        dimDensity,
        0.7
    );

    a_ = dimensionedScalar::lookupOrDefault
    (
        "a",
        BokilBewtraCoeffs_,
        dimless/Xc_.dimensions(),
        0.132
    );

    fac_ = BokilBewtraCoeffs_.lookupOrDefault<scalar>("Xscaling", 1.0);

    return true;
}


// ************************************************************************* //
