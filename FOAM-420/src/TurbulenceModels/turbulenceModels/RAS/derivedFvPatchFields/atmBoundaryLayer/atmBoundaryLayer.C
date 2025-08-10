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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "atmBoundaryLayer.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayer::atmBoundaryLayer()
:
    flowDir_(Zero),
    zDir_(Zero),
    kappa_(0.41),
    Cmu_(0.09),
    Uref_(0),
    Zref_(0),
    z0_(0),
    zGround_(0),
    Ustar_(0)
{}


atmBoundaryLayer::atmBoundaryLayer(const vectorField& p, const dictionary& dict)
:
    flowDir_(dict.lookup("flowDir")),
    zDir_(dict.lookup("zDir")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Uref_(readScalar(dict.lookup("Uref"))),
    Zref_(readScalar(dict.lookup("Zref"))),
    z0_("z0", dict, p.size()),
    zGround_("zGround", dict, p.size()),
    Ustar_(p.size())
{
    if (mag(flowDir_) < SMALL || mag(zDir_) < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    // Ensure direction vectors are normalized
    flowDir_ /= mag(flowDir_);
    zDir_ /= mag(zDir_);

    Ustar_ = kappa_*Uref_/(log((Zref_ + z0_)/z0_));
}


atmBoundaryLayer::atmBoundaryLayer
(
    const atmBoundaryLayer& ptf,
    const fvPatchFieldMapper& mapper
)
:
    flowDir_(ptf.flowDir_),
    zDir_(ptf.zDir_),
    kappa_(ptf.kappa_),
    Cmu_(ptf.Cmu_),
    Uref_(ptf.Uref_),
    Zref_(ptf.Zref_),
    z0_(mapper(ptf.z0_)),
    zGround_(mapper(ptf.zGround_)),
    Ustar_(mapper(ptf.Ustar_))
{}


atmBoundaryLayer::atmBoundaryLayer(const atmBoundaryLayer& blpvf)
:
    flowDir_(blpvf.flowDir_),
    zDir_(blpvf.zDir_),
    kappa_(blpvf.kappa_),
    Cmu_(blpvf.Cmu_),
    Uref_(blpvf.Uref_),
    Zref_(blpvf.Zref_),
    z0_(blpvf.z0_),
    zGround_(blpvf.zGround_),
    Ustar_(blpvf.Ustar_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayer::autoMap(const fvPatchFieldMapper& m)
{
    m(z0_, z0_);
    m(zGround_, zGround_);
    m(Ustar_, Ustar_);
}


void atmBoundaryLayer::rmap
(
    const atmBoundaryLayer& blptf,
    const labelList& addr
)
{
    z0_.rmap(blptf.z0_, addr);
    zGround_.rmap(blptf.zGround_, addr);
    Ustar_.rmap(blptf.Ustar_, addr);
}


void atmBoundaryLayer::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mapper.map(z0_, scalar(0));
    mapper.map(zGround_, scalar(0));
    mapper.map(Ustar_, scalar(0));
}


tmp<vectorField> atmBoundaryLayer::U(const vectorField& p) const
{
    scalarField Un
    (
        (Ustar_/kappa_)
       *log(((zDir_ & p) - zGround_ + z0_)/z0_)
    );

    return flowDir_*Un;
}


tmp<scalarField> atmBoundaryLayer::k(const vectorField& p) const
{
    return sqr(Ustar_)/sqrt(Cmu_);
}


tmp<scalarField> atmBoundaryLayer::epsilon(const vectorField& p) const
{
    return pow3(Ustar_)/(kappa_*((zDir_ & p) - zGround_ + z0_));
}


void atmBoundaryLayer::write(Ostream& os) const
{
    z0_.writeEntry("z0", os);
    os.writeEntry("flowDir", flowDir_);
    os.writeEntry("zDir", zDir_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("Uref", Uref_);
    os.writeEntry("Zref", Zref_);
    zGround_.writeEntry("zGround", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
