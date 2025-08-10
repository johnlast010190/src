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

\*---------------------------------------------------------------------------*/

#include "tetherPotential/derived/restrainedHarmonicSpring/restrainedHarmonicSpring.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(restrainedHarmonicSpring, 0);

addToRunTimeSelectionTable
(
    tetherPotential,
    restrainedHarmonicSpring,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

restrainedHarmonicSpring::restrainedHarmonicSpring
(
    const word& name,
    const dictionary& tetherPotentialProperties
)
:
    tetherPotential(name, tetherPotentialProperties),
    restrainedHarmonicSpringCoeffs_
    (
        tetherPotentialProperties.subDict(typeName + "Coeffs")
    ),
    springConstant_
    (
        readScalar(restrainedHarmonicSpringCoeffs_.lookup("springConstant"))
    ),
    rR_
    (
        readScalar(restrainedHarmonicSpringCoeffs_.lookup("rR"))
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar restrainedHarmonicSpring::energy(const vector r) const
{
    scalar magR = mag(r);

    if (magR < rR_)
    {
        return 0.5 * springConstant_ * magSqr(r);
    }
    else
    {
        return 0.5 * springConstant_ * rR_ * rR_
          + springConstant_ * rR_ * (magR - rR_);
    }
}


vector restrainedHarmonicSpring::force(const vector r) const
{
    scalar magR = mag(r);

    if (magR < rR_)
    {
        return -springConstant_ * r;
    }
    else
    {
        return -springConstant_ * rR_ * r / magR;
    }
}


bool restrainedHarmonicSpring::read
(
    const dictionary& tetherPotentialProperties
)
{
    tetherPotential::read(tetherPotentialProperties);

    restrainedHarmonicSpringCoeffs_ =
        tetherPotentialProperties.subDict(typeName + "Coeffs");

    restrainedHarmonicSpringCoeffs_.lookup("springConstant") >> springConstant_;
    restrainedHarmonicSpringCoeffs_.lookup("rR") >> rR_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// ************************************************************************* //
