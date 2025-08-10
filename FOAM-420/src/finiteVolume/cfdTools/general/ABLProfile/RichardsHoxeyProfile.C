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

#include "cfdTools/general/ABLProfile/RichardsHoxeyProfile.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace ABLProfiles
{

defineTypeNameAndDebug(RichardsHoxeyProfile, 0);
addToRunTimeSelectionTable(ABLProfile, RichardsHoxeyProfile, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RichardsHoxeyProfile::RichardsHoxeyProfile
(
    const fvPatch& patch
)
:
    ABLProfile(patch),
    Cmu_(0.0),
    Kappa_(0.0),
    Uref_(0.0),
    Href_(0.0),
    z0_(0.0),
    fLxyz_(vector::one)
{}


RichardsHoxeyProfile::RichardsHoxeyProfile
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
    Kappa_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<scalar>("Kappa", 0.41)
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
    z0_
    (
        readScalar
        (
            dict.subDict(this->typeName + "Coeffs").lookup("z0")
        )
    ),
    fLxyz_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<vector>("fLxyz", vector::one)
    )
{}


RichardsHoxeyProfile::RichardsHoxeyProfile
(
    const RichardsHoxeyProfile& abl
)
:
    ABLProfile(abl),
    Cmu_(abl.Cmu_),
    Kappa_(abl.Kappa_),
    Uref_(abl.Uref_),
    Href_(abl.Href_),
    z0_(abl.z0_),
    fLxyz_(abl.fLxyz_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void RichardsHoxeyProfile::setProfiles
(
    vectorField& U,
    scalarField& L,
    symmTensorField& R
)
{
    // Limit surfacxe roughness to 0.01 mm to prevent singularity
    const scalar z0min(max(z0_, 0.00001));

    const vector Ustar = Kappa_*Uref_/(log((Href_  + z0min)/z0min));

    U = (Ustar/Kappa_)*log((Z_->value() + z0min)/z0min);

    const scalarField k(U.size(), magSqr(Ustar)/sqrt(Cmu_));
    const scalarField epsilon(pow3(mag(Ustar))/(Kappa_*(Z_->value() + z0min)));
    const vectorField dUdz(Ustar*(1.0/(Kappa_*(Z_->value() + z0min))));
    const scalarField nut(Cmu_*sqr(k)/epsilon);
    R = 2.0/3.0*I*k - nut*symm(Z_->dir()*dUdz);

    L = k*sqrt(k)/epsilon*pow(Cmu_, 0.75);

}


void RichardsHoxeyProfile::setProfiles
(
    vectorField& U,
    vectorField& L,
    symmTensorField& R
)
{
    // Limit surfacxe roughness to 0.01 mm to prevent singularity
    const scalar z0min(max(z0_, 0.00001));

    const vector Ustar = Kappa_*Uref_/(log((Href_  + z0min)/z0min));

    U = (Ustar/Kappa_)*log((Z_->value() + z0min)/z0min);

    const scalarField k(U.size(), magSqr(Ustar)/sqrt(Cmu_));
    const scalarField epsilon(pow3(mag(Ustar))/(Kappa_*(Z_->value() + z0min)));
    const vectorField dUdz(Ustar*(1.0/(Kappa_*(Z_->value() + z0min))));
    const scalarField nut(Cmu_*sqr(k)/epsilon);
    R = 2.0/3.0*I*k - nut*symm(Z_->dir()*dUdz);

    L = fLxyz_*k*sqrt(k)/epsilon*pow(Cmu_, 0.75);
}


void RichardsHoxeyProfile::write(Ostream& os) const
{
    ABLProfile::write(os);

    os.writeEntry("profileType", typeName);

    word profileCoeffs = typeName + word("Coeffs");
    os.beginBlock(profileCoeffs);
    if (Cmu_ != 0.09)
    {
        os.writeEntry("Cmu", Cmu_);
    }
    if (Kappa_ != 0.41)
    {
        os.writeEntry("Kappa", Kappa_);
    }
    os.writeEntry("Uref", Uref_);
    os.writeEntry("Href", Href_);
    os.writeEntry("z0", z0_);

    if (fLxyz_ != vector::one)
    {
        os.writeEntry("fLxyz", fLxyz_);
    }

    os.endBlock() << endl;
}

} // End namespace ABLProfiles
} // End namespace Foam

// ************************************************************************* //
