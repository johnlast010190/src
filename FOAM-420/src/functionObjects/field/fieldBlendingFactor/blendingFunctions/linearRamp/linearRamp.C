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
#include "fieldBlendingFactor/blendingFunctions/linearRamp/linearRamp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingFunctions
{

defineTypeNameAndDebug(linearRamp, 0);
addToRunTimeSelectionTable(blendingFunction, linearRamp, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearRamp::linearRamp(const dictionary& dict, bool dum)
:
    blendingFunction(dict),
    minValue_(readScalar(dict.lookup("minValue"))),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    minCoeff_(readScalar(dict.lookup("minCoeff"))),
    maxCoeff_(readScalar(dict.lookup("maxCoeff"))),
    linearF_()
{
    linearF_[1] = (minCoeff_ - maxCoeff_)/(minValue_ - maxValue_);
    linearF_[0] = minCoeff_ - linearF_[1]*minValue_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void linearRamp::update(scalar value)
{
    minValue_ = value;
    linearF_[1] = (minCoeff_ - maxCoeff_)/(minValue_ - maxValue_);
    linearF_[0] = minCoeff_ - linearF_[1]*minValue_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

scalar linearRamp::operator()(scalar value) const
{
    if (value > maxValue_)
    {
        return maxCoeff_;
    }
    else if (value < minValue_)
    {
        return minCoeff_;
    }

    return linearF_.value(value);
}


} // End namespace blendingFunctions
} // End namespace Foam

// ************************************************************************* //
