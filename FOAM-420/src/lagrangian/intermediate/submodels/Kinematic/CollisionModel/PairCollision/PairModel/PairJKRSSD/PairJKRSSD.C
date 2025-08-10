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

#include "submodels/Kinematic/CollisionModel/PairCollision/PairModel/PairJKRSSD/PairJKRSSD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::PairJKRSSD<CloudType>::findMinMaxProperties
(
    scalar& RMin,
    scalar& rhoMax,
    scalar& UMagMax
) const
{
    RMin = VGREAT;
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

        RMin = min(dEff, RMin);

        rhoMax = max(p.rho(), rhoMax);

        UMagMax = max
        (
            mag(p.U()) + mag(p.omega())*dEff/2,
            UMagMax
        );
    }

    // Transform the minimum diameter into minimum radius
    //     rMin = dMin/2
    // then rMin into minimum R,
    //     1/RMin = 1/rMin + 1/rMin,
    //     RMin = rMin/2 = dMin/4

    RMin /= 4.0;

    // Multiply by two to create the worst-case relative velocity

    UMagMax = 2*UMagMax;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairJKRSSD<CloudType>::PairJKRSSD
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PairModel<CloudType>(dict, cloud, typeName),
    Estar_(),
    Gstar_(),
    alpha_(readScalar(this->coeffDict().lookup("alpha"))),
    mu_(readScalar(this->coeffDict().lookup("mu"))),
    gamma_(readScalar(this->coeffDict().lookup("gamma"))),
    collisionResolutionSteps_
    (
        readScalar(this->coeffDict().lookup("collisionResolutionSteps"))
    ),
    volumeFactor_(1.0),
    useEquivalentSize_(Switch(this->coeffDict().lookup("useEquivalentSize"))),
    verbose_(Switch(this->coeffDict().template lookupOrDefault("verbose",false)))
{
    if (useEquivalentSize_)
    {
        volumeFactor_ = readScalar(this->coeffDict().lookup("volumeFactor"));
    }

    scalar nu = this->owner().constProps().poissonsRatio();

    scalar E = this->owner().constProps().youngsModulus();

    Estar_ = E/(2.0*(1.0 - sqr(nu)));

    scalar G = E/(2.0*(1.0 + nu));

    Gstar_ = G/(2.0*(2.0 - nu));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairJKRSSD<CloudType>::~PairJKRSSD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PairJKRSSD<CloudType>::controlsTimestep() const
{
    return true;
}


template<class CloudType>
Foam::label Foam::PairJKRSSD<CloudType>::nSubCycles() const
{
    if (!(this->owner().size()))
    {
        return 1;
    }

    scalar RMin;
    scalar rhoMax;
    scalar UMagMax;

    findMinMaxProperties(RMin, rhoMax, UMagMax);

    // Note:  pi^(7/5)*(5/4)^(2/5) = 5.429675
    scalar minCollisionDeltaT =
        5.429675
       *RMin
       *pow(rhoMax/(Estar_*sqrt(UMagMax) + VSMALL), 0.4)
       /collisionResolutionSteps_;

    return ceil(this->owner().time().deltaTValue()/minCollisionDeltaT);
}


template<class CloudType>
void Foam::PairJKRSSD<CloudType>::evaluatePair
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    vector r_AB = (pA.position() - pB.position());

    scalar dAEff = pA.d();

    if (useEquivalentSize_)
    {
        dAEff *= cbrt(pA.nParticle()*volumeFactor_);
    }

    scalar dBEff = pB.d();

    if (useEquivalentSize_)
    {
        dBEff *= cbrt(pB.nParticle()*volumeFactor_);
    }

    scalar r_ABmag = mag(r_AB);

    scalar penDepth = 0.5*(dAEff + dBEff) - r_ABmag;

    if (penDepth > 0)
    {
        //Particles in collision

        vector rHat_AB = r_AB/(r_ABmag + VSMALL);

        vector U_AB = pA.U() - pB.U();

        // Effective radius
        scalar R = 0.5*dAEff*dBEff/(dAEff + dBEff);

        // Effective mass
        scalar M = pA.mass()*pB.mass()/(pA.mass() + pB.mass());

        scalar kN = (4.0/3.0)*sqrt(R)*Estar_;

        scalar etaN = alpha_*sqrt(M*kN)*pow025(penDepth);

        // Normal force
        vector fN_AB =
            rHat_AB
           *(kN*pow(penDepth, 3.0/2.0) - etaN*(U_AB & rHat_AB));

        // If the work of adhesion is non-zero
        if (gamma_)
        {
            fN_AB +=
                - 3.0*gamma_
                * constant::mathematical::pi*R
                * rHat_AB;
        }

        pA.f() += fN_AB;
        pB.f() += -fN_AB;
        pA.parcelForces() += fN_AB;
        pB.parcelForces() += -fN_AB;

        if (verbose_)
        {
            Info<<"     Parcel-Parcel Collision "<<endl;
            Info<<"     Effective radius "<<R<<endl;
            Info<<"     Penetration depth "<<penDepth<<endl;
            Info<<"     Normal Force "<<fN_AB<<endl;
        }

        vector USlip_AB =
            U_AB - (U_AB & rHat_AB)*rHat_AB
          + (pA.omega() ^ (dAEff/2*-rHat_AB))
          - (pB.omega() ^ (dBEff/2*rHat_AB));

        scalar deltaT = this->owner().mesh().time().deltaTValue();

        vector& tangentialOverlap_AB =
            pA.collisionRecords().matchPairRecord
            (
                pB.origProc(),
                pB.origId()
            ).collisionData();

        vector& tangentialOverlap_BA =
            pB.collisionRecords().matchPairRecord
            (
                pA.origProc(),
                pA.origId()
            ).collisionData();

        vector deltaTangentialOverlap_AB = USlip_AB*deltaT;

        tangentialOverlap_AB += deltaTangentialOverlap_AB;
        tangentialOverlap_BA += -deltaTangentialOverlap_AB;

        scalar tangentialOverlapMag = mag(tangentialOverlap_AB);

        if (tangentialOverlapMag > VSMALL)
        {
            scalar kT = 8.0*sqrt(R*penDepth)*Gstar_;

            scalar etaT = etaN;

            // Tangential force
            vector fT_AB;

            if (kT*tangentialOverlapMag > mu_*mag(fN_AB))
            {
                // Tangential force greater than sliding friction,
                // particle slips

                fT_AB = -mu_*mag(fN_AB)*USlip_AB/mag(USlip_AB);

                tangentialOverlap_AB = vector::zero;
                tangentialOverlap_BA = vector::zero;
            }
            else
            {
                fT_AB =
                    -kT*tangentialOverlapMag
                   *tangentialOverlap_AB/tangentialOverlapMag
                  - etaT*USlip_AB;
            }

            pA.f() += fT_AB;
            pB.f() += -fT_AB;
            pA.parcelForces() += fT_AB;
            pB.parcelForces() += -fT_AB;

            pA.torque() += (dAEff/2*-rHat_AB) ^ fT_AB;
            pB.torque() += (dBEff/2*rHat_AB) ^ -fT_AB;
        }
    }
}


// ************************************************************************* //
