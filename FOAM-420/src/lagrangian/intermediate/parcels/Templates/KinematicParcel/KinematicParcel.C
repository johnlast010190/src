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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "parcels/Templates/KinematicParcel/KinematicParcel.H"
#include "submodels/Kinematic/ParticleForces/forceSuSp/forceSuSp.H"
#include "IntegrationScheme/IntegrationScheme/IntegrationScheme.H"
#include "meshTools/meshTools.H"
#include "clouds/Templates/KinematicCloud/cloudSolution/cloudSolution.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::KinematicParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    tetIndices tetIs = this->currentTetIndices();

    rhoc_ = td.rhoInterp().interpolate(this->position(), tetIs);

    if (rhoc_ < td.cloud().constProps().rhoMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed density in cell " << celli << " to "
                << td.cloud().constProps().rhoMin() <<  nl << endl;
        }

        rhoc_ = td.cloud().constProps().rhoMin();
    }

    Uc_ = td.UInterp().interpolate(this->position(), tetIs);

    muc_ = td.muInterp().interpolate(this->position(), tetIs);

    // Apply dispersion components to carrier phase velocity
    Uc_ = td.cloud().dispersion().update
    (
        dt,
        celli,
        U_,
        Uc_,
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    Uc_ += td.cloud().UTrans()[celli]/massCell(celli);
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    scalar Re = this->Re(U_, d_, rhoc_, muc_);

    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Motion
    // ~~~~~~

    // instantaneous parcel forces
    vector pForces = vector::zero;

    // Calculate new particle velocity
    this->U_ = calcVelocity(td, dt, celli, Re, muc_, mass0, Su, dUTrans, Spu, pForces);

    // update the parcel forces
    this->parcelForces_ = pForces;

    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().solution().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[celli] += np0*dUTrans;

        // Update momentum transfer coefficient
        td.cloud().UCoeff()[celli] += np0*Spu;
    }
}


template<class ParcelType>
template<class TrackData>
const Foam::vector Foam::KinematicParcel<ParcelType>::calcVelocity
(
    TrackData& td,
    const scalar dt,
    const label celli,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu,
    vector& pForces
) const
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;

    const forceType& forces = td.cloud().forces();

    // Momentum source due to particle forces
    const parcelType& p = static_cast<const parcelType&>(*this);
    const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, mu);
    const forceSuSp Fncp = forces.calcNonCoupled(p, dt, mass, Re, mu);
    const forceSuSp Feff = Fcp + Fncp;
    const scalar massEff = forces.massEff(p, mass);


    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const vector abp = (Feff.Sp()*Uc_ + (Feff.Su() + Su))/massEff;
    const scalar bp = Feff.Sp()/massEff;

    // keep a record of the forces on this current parcel
    pForces = (Feff.Sp()*Uc_ + (Feff.Su() + Su));

    Spu = dt*Feff.Sp();

    IntegrationScheme<vector>::integrationResult Ures =
        td.cloud().UIntegrator().integrate(U_, dt, abp, bp);

    vector Unew = Ures.value();

    // note: Feff.Sp() and Fc.Sp() must be the same
    dUTrans += dt*(Feff.Sp()*(Ures.average() - Uc_) - Fcp.Su());

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = td.cloud().pMesh();
    meshTools::constrainDirection(mesh, mesh.geometricD(), Unew);
    meshTools::constrainDirection(mesh, mesh.geometricD(), dUTrans);

    // Apply restriction if cell is smaller than parcel
    scalar cellVolume = mesh.cellVolumes()[celli];
    scalar parcelVolume = p.volume();

    if (
        td.cloud().solution().decoupleSmallCells()
        && (parcelVolume > cellVolume)
    )
    {
        if (debug)
        {
            Pout
                << "    Parcel volume " << parcelVolume << " is larger than "
                << " cell "<< celli <<" volume of "<< cellVolume << nl
                << "     decoupling momentum sources in this cell."<< endl;
        }

        dUTrans = vector(0,0,0);
        Spu = scalar(0);
    }


    return Unew;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p
)
:
    ParcelType(p),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    parcelID_(p.parcelID_),
    parcelForces_(p.parcelForces_),
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_)
{}


template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    parcelID_(p.parcelID_),
    parcelForces_(p.parcelForces_),
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
    const cloudSolution& solution = td.cloud().solution();
    const scalarField& cellLengthScale = td.cloud().cellLengthScale();

    scalar tEnd = (1.0 - p.stepFraction())*trackTime;
    scalar dtMax = solution.deltaTMax(trackTime);

    bool tracking = true;
    label nTrackingStalled = 0;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Apply correction to position for reduced-D cases
        meshTools::constrainToMeshCentre(mesh, p.position());

        const point start(p.position());

        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Cache the parcel current cell as this will change if a face is hit
        const label celli = p.cell();

        const scalar magU = mag(U_);
        if (p.active() && tracking && (magU > ROOTVSMALL))
        {
            const scalar d = dt*magU;
            const scalar deltaLMax = solution.deltaLMax(cellLengthScale[celli]);
            const scalar dCorr = min(d, deltaLMax);
            dt *=
                dCorr/d
               *p.trackToFace(p.position() + dCorr*U_/magU, td);
        }

        tEnd -= dt;

        const scalar newStepFraction = 1.0 - tEnd/trackTime;

        if (tracking)
        {
            if
            (
                mag(p.stepFraction() - newStepFraction)
              < particle::minStepFractionTol
            )
            {
                nTrackingStalled++;

                if (nTrackingStalled > maxTrackAttempts)
                {
                    tracking = false;
                }
            }
            else
            {
                nTrackingStalled = 0;
            }
        }

        p.stepFraction() = newStepFraction;

        bool calcParcel = true;
        if (!tracking && solution.steadyState())
        {
            calcParcel = false;
        }

        // Avoid problems with extremely small timesteps
        if ((dt > ROOTVSMALL) && calcParcel)
        {
            // Update cell based properties
            p.setCellValues(td, dt, celli);

            if (solution.cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(td, dt, celli);
            }

            p.calc(td, dt, celli);
        }

        if (p.onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
            {
                td.switchProcessor = true;
            }
        }

        p.age() += dt;

        td.cloud().functions().postMove(p, celli, dt, start, td.keepParticle);
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitFace(TrackData& td)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.cloud().functions().postFace(p, p.face(), td.keepParticle);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitFace(int& td)
{}


template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch& pp,
    TrackData& td,
    const label patchi,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    // Invoke post-processing model
    td.cloud().functions().postPatch
    (
        p,
        pp,
        trackFraction,
        tetIs,
        td.keepParticle
    );

    if (isA<processorPolyPatch>(pp))
    {
        // Skip processor patches
        return false;
    }
    else if (td.cloud().surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        // Surface film model consumes the interaction, i.e. all
        // interactions done
        return true;
    }
    else
    {
        // Invoke patch interaction model
        return td.cloud().patchInteraction().correct
        (
            p,
            pp,
            td.keepParticle,
            trackFraction,
            tetIs
        );
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td,
    const tetIndices&
)
{
    // Wall interactions handled by generic hitPatch function
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;

    td.cloud().patchInteraction().addToEscapedParcels(nParticle_*mass());
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const transformer& transform
)
{
    ParcelType::transformProperties(transform);

    U_ = transform.transform(U_);
}


template<class ParcelType>
Foam::scalar Foam::KinematicParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    return 0.5*d_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "parcels/Templates/KinematicParcel/KinematicParcelIO.C"

// ************************************************************************* //
