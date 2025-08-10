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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/ABLProfile/AIJProfile.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace ABLProfiles
{

defineTypeNameAndDebug(AIJProfile, 0);
addToRunTimeSelectionTable(ABLProfile, AIJProfile, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AIJProfile::AIJProfile
(
    const fvPatch& patch
)
:
    ABLProfile(patch),
    Cmu_(0.0),
    Uref_(vector::zero),
    Href_(0.0),
    alpha_(0.0),
    zG_(0.0),
    AIJPaper_(false),
    fLxyz_(vector::one)
{}


AIJProfile::AIJProfile
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    ABLProfile(patch, dict),
    Cmu_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<scalar>("Cmu", 0.09)
    ),
    Uref_
    (
        dict.subDict(this->typeName + "Coeffs").lookup("Uref")
    ),
    Href_
    (
        readScalar
        (
            dict.subDict(this->typeName + "Coeffs").lookup("Href")
        )
    ),
    alpha_
    (
        readScalar
        (
            dict.subDict(this->typeName + "Coeffs").lookup("alpha")
        )
    ),
    zG_
    (
        readScalar
        (
            dict.subDict(this->typeName + "Coeffs").lookup("zG")
        )
     ),
    AIJPaper_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<bool>("consistentTKE", false)
    ),
    fLxyz_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<vector>("fLxyz", vector::one)
    )
{}

AIJProfile::AIJProfile
(
    const AIJProfile& abl
)
:
    ABLProfile(abl),
    Cmu_(abl.Cmu_),
    Uref_(abl.Uref_),
    Href_(abl.Href_),
    alpha_(abl.alpha_),
    zG_(abl.zG_),
    AIJPaper_(abl.AIJPaper_),
    fLxyz_(abl.fLxyz_)
{}

// * * * * * * * * * * * * * * * * * Functions* * * * * * * * * * * * * * * //

void AIJProfile::setProfiles
(
    vectorField& U,
    scalarField& L,
    symmTensorField& R
)
{
    U = Uref_*pow(Z_->value()/Href_, alpha_);

    scalar fTKE;
    if (AIJPaper_)
    {
        fTKE = 1.0;
    }
    else
    {
        fTKE = 3.0/2.0;
    }

    scalarField k(fTKE*sqr(0.1*pow(Z_->value()/zG_, -alpha_-0.05)*mag(U)));

    scalarField epsilon
    (
        sqrt(Cmu_)*k*mag(Uref_)/Href_*alpha_*pow(Z_->value()/Href_, alpha_-1.0)
    );

    scalarField nut ( Cmu_*sqr(k)/epsilon );

    vectorField dUdz(alpha_*Uref_*pow(Z_->value()/Href_, alpha_)/Z_->value());

    L = k*sqrt(k)/epsilon * pow(Cmu_,0.75);

    R = 2.0/3.0*I*k - nut*symm(Z_->dir() * dUdz);

}

void AIJProfile::setProfiles
(
    vectorField& U,
    vectorField& L,
    symmTensorField& R
)
{
    U = Uref_*pow(Z_->value()/Href_, alpha_);

    scalar fTKE;
    if (AIJPaper_)
    {
        fTKE = 1.0;
    }
    else
    {
        fTKE = 3.0/2.0;
    }

    scalarField k(fTKE*sqr(0.1*pow(Z_->value()/zG_, -alpha_-0.05)*mag(U)));

    scalarField epsilon
    (
        sqrt(Cmu_)*k*mag(Uref_)/Href_*alpha_*pow(Z_->value()/Href_, alpha_-1.0)
    );

    scalarField nut ( Cmu_*sqr(k)/epsilon );

    vectorField dUdz(alpha_*Uref_*pow(Z_->value()/Href_, alpha_)/Z_->value());

    L = fLxyz_*k*sqrt(k)/epsilon * pow(Cmu_,0.75);

    R = 2.0/3.0*I*k - nut*symm(Z_->dir() * dUdz);

}


void AIJProfile::write(Ostream& os) const
{
    ABLProfile::write(os);

    os.writeEntry("profileType", typeName);

    word profileCoeffs = typeName + word("Coeffs");
    os.beginBlock(profileCoeffs);
    if (Cmu_ != 0.09)
    {
        os.writeEntry("Cmu", Cmu_);
    }
    os.writeEntry("Uref", Uref_);
    os.writeEntry("Href", Href_);
    os.writeEntry("alpha", alpha_);
    os.writeEntry("zG", zG_);
    os.writeEntry("consistentTKE", AIJPaper_);

    if (fLxyz_!=vector::one)
    {
        os.writeEntry("fLxyz", fLxyz_);
    }

    os.endBlock() << endl;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //



} // End namespace ABLProfiles
} // End namespace Foam

// ************************************************************************* //
