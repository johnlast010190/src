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
    (c) 2011 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingFunctions/heaviside/heaviside.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingFunctions
{

defineTypeNameAndDebug(heaviside, 0);
addToRunTimeSelectionTable(blendingFunction, heaviside, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

heaviside::heaviside(const dictionary& dict, bool dum)
:
    blendingFunction(dict),
    posCoeff_(readScalar(dict.lookup("posCoeff"))),
    negCoeff_(readScalar(dict.lookup("negCoeff"))),
    intercept_(dict.lookupOrDefault<scalar>("interceptValue", 0))
{}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

scalar heaviside::operator()(scalar value) const
{
    value -= intercept_;

    if (value > 0)
    {
        return posCoeff_;
    }
    else if (value < 0)
    {
        return negCoeff_;
    }

    return 0.5*(posCoeff_ + negCoeff_);
}


} // End namespace blendingFunctions
} // End namespace Foam

// ************************************************************************* //
