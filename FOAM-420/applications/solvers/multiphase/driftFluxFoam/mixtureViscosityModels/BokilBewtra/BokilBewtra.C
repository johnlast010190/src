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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "BokilBewtra.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(BokilBewtra, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        BokilBewtra,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::BokilBewtra::BokilBewtra
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityProperties, U, phi),
    BokilBewtraCoeffs_(viscosityProperties.optionalSubDict(modelName + "Coeffs")),
    BokilBewtraViscosityCoeff_
    (
        "coeff",
        dimensionSet(1, -1, -1, 0, 0),
        BokilBewtraCoeffs_.lookup("coeff")
    ),
    BokilBewtraViscosityExponent_
    (
        "exponent",
        dimless,
        BokilBewtraCoeffs_.lookup("exponent")
    ),
    muMax_
    (
        "muMax",
        dimensionSet(1, -1, -1, 0, 0),
        BokilBewtraCoeffs_.lookup("muMax")
    ),
    alphac_
    (
        "alphac",
        dimensionSet(0, 0, 0, 0, 0),
        BokilBewtraCoeffs_.lookup("alphac")
    ),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityProperties.lookupOrDefault<word>("alpha", "alpha"),
                viscosityProperties.dictName()
            )
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::BokilBewtra::mu(const volScalarField& muc) const
{
    volScalarField blending
        (
            min
                (
                    max((alphac_ - alpha_) / alphac_, scalar(0.0)), scalar(1.0)
                )
        );

    return min
    (
        muc * blending
      + BokilBewtraViscosityCoeff_ * (1.-blending)
        * pow
        (
            scalar(10),
            BokilBewtraViscosityExponent_*alpha_
        ),
        muMax_
    );
}


bool Foam::mixtureViscosityModels::BokilBewtra::read
(
    const dictionary& viscosityProperties
)
{
    mixtureViscosityModel::read(viscosityProperties);

    BokilBewtraCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    BokilBewtraCoeffs_.lookup("coeff") >> BokilBewtraViscosityCoeff_;
    BokilBewtraCoeffs_.lookup("exponent") >> BokilBewtraViscosityExponent_;
    BokilBewtraCoeffs_.lookup("muMax") >> muMax_;
    BokilBewtraCoeffs_.lookup("alphac") >> alphac_;

    return true;
}


// ************************************************************************* //
