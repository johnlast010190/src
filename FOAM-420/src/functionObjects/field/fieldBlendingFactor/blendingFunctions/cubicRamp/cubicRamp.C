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
#include "fieldBlendingFactor/blendingFunctions/cubicRamp/cubicRamp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingFunctions
{

defineTypeNameAndDebug(cubicRamp, 0);
addToRunTimeSelectionTable(blendingFunction, cubicRamp, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cubicRamp::cubicRamp(const dictionary& dict, bool dum)
:
    blendingFunction(dict),
    minValue_(readScalar(dict.lookup("minValue"))),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    minCoeff_(readScalar(dict.lookup("minCoeff"))),
    maxCoeff_(readScalar(dict.lookup("maxCoeff"))),
    cubicF_()
{

    scalar x1 = minValue_;
    scalar x12 = sqr(x1);
    scalar x13 = x12*x1;
    scalar x2 = maxValue_;
    scalar x22 = sqr(x2);
    scalar x23 = x22*x2;

    scalar cf3Rcp = 3*x12*x2 - x23 -2*x13
        - 3.0/2.0 *(x12-x22)*(2*x1*x2-x22-x12)/(x1 - x2);

    cubicF_[3] = (minCoeff_ - maxCoeff_) / cf3Rcp;
        //= -(minCoeff_ - maxCoeff_)/(2*(::pow(minValue_,3) - ::pow(maxValue_,3))
        //- 3/2*(minValue_ + maxValue_)*(sqr(minValue_)-sqr(maxValue_)));

    cubicF_[2] = -3.0/2.0*cubicF_[3]*(x12-x22)/(x1-x2);
        //-3/2*cubicF_[3]*(minValue_ + maxValue_);

    cubicF_[1] = -2*cubicF_[2]*x2 - 3*cubicF_[3]*x22;
        //-3*cubicF_[3]*sqr(minValue_)
        //- 2*cubicF_[2]*minValue_;

    cubicF_[0] = minCoeff_ - cubicF_[1]*x1 - cubicF_[2]*x12 - cubicF_[3]*x13;
        //minCoeff_ - cubicF_[3]*::pow(minValue_,3)
        //- cubicF_[2]*sqr(minValue_) - cubicF_[1]*minValue_;

    if (debug)
    {
        Info<< "Input : minValue - " << minValue_ << ", maxValue - "
             << maxValue_ << ", minCoeff - " << minCoeff_
             << ", maxCoeff - " << maxCoeff_ << endl;

        Info<< "Cubic ramp coefficients: " << cubicF_ << endl;
        Info<< "Return min: " << this->operator()(minValue_) << endl;
        Info<< "Return mean: " << this->operator()(0.5*(minValue_+maxValue_)) << endl;
        Info<< "Return max: " << this->operator()(maxValue_) << endl;
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void cubicRamp::update(scalar value)
{
    minValue_ = value;

    scalar x1 = minValue_;
    scalar x12 = sqr(x1);
    scalar x13 = x12*x1;
    scalar x2 = maxValue_;
    scalar x22 = sqr(x2);
    scalar x23 = x22*x2;

    scalar cf3Rcp = 3*x12*x2 - x23 -2*x13
        - 3.0/2.0 *(x12-x22)*(2*x1*x2-x22-x12)/(x1 - x2);
    cubicF_[3] = (minCoeff_ - maxCoeff_) / cf3Rcp;
    cubicF_[2] = -3.0/2.0*cubicF_[3]*(x12-x22)/(x1-x2);
    cubicF_[1] = -2*cubicF_[2]*x2 - 3*cubicF_[3]*x22;
    cubicF_[0] = minCoeff_ - cubicF_[1]*x1 - cubicF_[2]*x12 - cubicF_[3]*x13;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

scalar cubicRamp::operator()(scalar value) const
{
    if (value > maxValue_)
    {
        return maxCoeff_;
    }
    else if (value < minValue_)
    {
        return minCoeff_;
    }

    return cubicF_.value(value);
}


} // End namespace blendingFunctions
} // End namespace Foam

// ************************************************************************* //
