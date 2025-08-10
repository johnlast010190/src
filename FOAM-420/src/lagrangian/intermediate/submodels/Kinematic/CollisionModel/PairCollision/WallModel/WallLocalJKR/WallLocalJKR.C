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

#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallLocalJKR/WallLocalJKR.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::WallLocalJKR<CloudType>::findMinMaxProperties
(
    scalar& rMin,
    scalar& rhoMax,
    scalar& UMagMax
) const
{
    rMin = VGREAT;
    rhoMax = -VGREAT;
    UMagMax = -VGREAT;

    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        const typename CloudType::parcelType& p = iter();

        // Finding minimum diameter to avoid excessive arithmetic

        scalar dEff = p.d();

        if (useEquivalentSize_)
        {
            dEff *= cbrt(p.nParticle()*volumeFactor_);
        }

        rMin = min(dEff, rMin);

        rhoMax = max(p.rho(), rhoMax);

        UMagMax = max
        (
            mag(p.U()) + mag(p.omega())*dEff/2,
            UMagMax
        );
    }

    // Transform the minimum diameter into minimum radius
    //     rMin = dMin/2

    rMin /= 2.0;
}


template<class CloudType>
void Foam::WallLocalJKR<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const point& site,
    const WallSiteData<vector>& data,
    scalar pREff,
    bool cohesion
) const
{
    // wall patch index
    label wPI = patchMap_[data.patchIndex()];

    // data for this patch
    scalar Estar = Estar_[wPI];

    scalar gamma = gamma_[wPI];

    vector r_PW = p.position() - site;

    scalar r_PW_mag = mag(r_PW);

    scalar penDepth = max(pREff - r_PW_mag, 0.0);

    vector rHat_PW = r_PW/(r_PW_mag + VSMALL);

    scalar kN = (4.0/3.0)*sqrt(pREff)*Estar;

    vector fN_PW =
        rHat_PW
       *(kN*pow(penDepth, 3.0/2.0));

    // If the work of adhesion is non-zero
    if (gamma > 0 && (penDepth>=0))
    {
        fN_PW +=
                -rHat_PW
               * (
                    3.0/2.0*gamma*constant::mathematical::pi*pREff
                 );
    }

    p.f() += fN_PW;
    p.parcelForces() += fN_PW;

    if (verbose_)
    {
        Info<<"     Parcel-Patch Collision "<<endl;
        Info<<"     Effective radius "<<pREff<<endl;
        Info<<"     Penetration depth "<<penDepth<<endl;
        Info<<"     Normal Force "<<fN_PW<<endl;
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WallLocalJKR<CloudType>::WallLocalJKR
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallModel<CloudType>(dict, cloud, typeName),
    Estar_(),
    gamma_(),
    patchMap_(),
    maxEstarIndex_(-1),
    collisionResolutionSteps_
    (
        readScalar
        (
            this->coeffDict().lookup("collisionResolutionSteps")
        )
    ),
    volumeFactor_(1.0),
    useEquivalentSize_(Switch(this->coeffDict().lookup("useEquivalentSize"))),
    verbose_(Switch(this->coeffDict().template lookupOrDefault("verbose",false)))
{
    if (useEquivalentSize_)
    {
        volumeFactor_ = readScalar(this->coeffDict().lookup("volumeFactor"));
    }

    scalar pNu = this->owner().constProps().poissonsRatio();

    scalar pE = this->owner().constProps().youngsModulus();

    const polyMesh& mesh = cloud.mesh();

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    patchMap_.setSize(bMesh.size(), -1);

    DynamicList<label> wallPatchIndices;

    forAll(bMesh, patchI)
    {
        if (isA<wallPolyPatch>(bMesh[patchI]))
        {
            wallPatchIndices.append(bMesh[patchI].index());
        }
    }

    label nWallPatches = wallPatchIndices.size();

    Estar_.setSize(nWallPatches);
    gamma_.setSize(nWallPatches);

    scalar maxEstar = -GREAT;

    forAll(wallPatchIndices, wPI)
    {
        const dictionary& patchCoeffDict
        (
            this->coeffDict().subDict(bMesh[wallPatchIndices[wPI]].name())
        );

        patchMap_[wallPatchIndices[wPI]] = wPI;

        scalar nu = readScalar(patchCoeffDict.lookup("poissonsRatio"));

        scalar E = readScalar(patchCoeffDict.lookup("youngsModulus"));

        scalar patchGamma = readScalar(patchCoeffDict.lookup("patchGamma"));

        scalar parcelGamma = readScalar(patchCoeffDict.lookup("parcelGamma"));

        gamma_[wPI] = patchGamma + parcelGamma;

        Estar_[wPI] = 1/((1 - sqr(pNu))/pE + (1 - sqr(nu))/E);

        if (Estar_[wPI] > maxEstar)
        {
            maxEstarIndex_ = wPI;

            maxEstar = Estar_[wPI];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WallLocalJKR<CloudType>::~WallLocalJKR()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::WallLocalJKR<CloudType>::pREff
(
    const typename CloudType::parcelType& p
) const
{
    if (useEquivalentSize_)
    {
        return p.d()/2*cbrt(p.nParticle()*volumeFactor_);
    }
    else
    {
        return p.d()/2;
    }
}


template<class CloudType>
bool Foam::WallLocalJKR<CloudType>::controlsTimestep() const
{
    return true;
}


template<class CloudType>
Foam::label Foam::WallLocalJKR<CloudType>::nSubCycles() const
{
    if (!(this->owner().size()))
    {
        return 1;
    }

    scalar rMin;
    scalar rhoMax;
    scalar UMagMax;

    findMinMaxProperties(rMin, rhoMax, UMagMax);

    // Note:  pi^(7/5)*(5/4)^(2/5) = 5.429675
    scalar minCollisionDeltaT =
        5.429675
       *rMin
       *pow(rhoMax/(Estar_[maxEstarIndex_]*sqrt(UMagMax) + VSMALL), 0.4)
       /collisionResolutionSteps_;

    return ceil(this->owner().time().deltaTValue()/minCollisionDeltaT);
}


template<class CloudType>
void Foam::WallLocalJKR<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSitePoints,
    const List<WallSiteData<vector>>& flatSiteData,
    const List<point>& sharpSitePoints,
    const List<WallSiteData<vector>>& sharpSiteData
) const
{
    scalar pREff = this->pREff(p);

    forAll(flatSitePoints, siteI)
    {
        evaluateWall
        (
            p,
            flatSitePoints[siteI],
            flatSiteData[siteI],
            pREff,
            true
        );
    }

    forAll(sharpSitePoints, siteI)
    {
        // Treating sharp sites like flat sites, except suppress cohesion

        evaluateWall
        (
            p,
            sharpSitePoints[siteI],
            sharpSiteData[siteI],
            pREff,
            false
        );
    }
}


// ************************************************************************* //
