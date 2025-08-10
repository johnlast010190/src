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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "panaRotorDiskSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "sources/derived/rotorDiskSource/trimModel/trimModel/trimModel.H"
#include "global/unitConversion/unitConversion.H"
#include "fvMatrices/fvMatrices.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(panaRotorDiskSource, 0);
        addToRunTimeSelectionTable(option, panaRotorDiskSource, dictionary);
    }

    template<> const char* NamedEnum<fv::panaRotorDiskSource::geometryModeType, 2>::
        names[] =
    {
        "auto",
        "specified"
    };

    const NamedEnum<fv::panaRotorDiskSource::geometryModeType, 2>
        fv::panaRotorDiskSource::geometryModeTypeNames_;

    template<> const char* NamedEnum<fv::panaRotorDiskSource::inletFlowType, 3>::
        names[] =
    {
        "fixed",
        "surfaceNormal",
        "local"
    };

    const NamedEnum<fv::panaRotorDiskSource::inletFlowType, 3>
        fv::panaRotorDiskSource::inletFlowTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::panaRotorDiskSource::checkData()
{
    // Set inflow type
    switch (selectionMode())
    {
        case smCellSet:
        case smCellZone:
        case smAll:
        {
            // Set the profile ID for each blade section
            profiles_.connectBlades(blade_.profileName(), blade_.profileID());
            switch (inletFlow_)
            {
                case ifFixed:
                {
                    coeffs_.lookup("inletVelocity") >> inletVelocity_;
                    break;
                }
                case ifSurfaceNormal:
                {
                    flowRate0_ =   readScalar(coeffs_.lookup("inletNormalVolumeFlux"));
                    // Convert Volume flux [ m^3/min ] to m/s
                    // In Pana system enforce flow along axis.
                    flowRate0_  /= 60.; // m^3/s
                    flowRate_ = flowRate0_;
                    fanAxis_ = coeffs_.lookup("axis");
                    fanAxis_ /= mag(fanAxis_);
                    inletVelocity_ = fanAxis_ * flowRate0_ / rotorDiskA_;
                    Info<<  "    - inletVelocity        = "<< inletVelocity_ << endl;

                    break;
                }
                case ifLocal:
                {
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown inlet velocity type" << abort(FatalError);
                }
            }


            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Source cannot be used with '"
                << selectionModeTypeNames_[selectionMode()]
                << "' mode.  Please use one of: " << nl
                << selectionModeTypeNames_[smCellSet] << nl
                << selectionModeTypeNames_[smCellZone] << nl
                << selectionModeTypeNames_[smAll]
                << exit(FatalError);
        }
    }
}


void Foam::fv::panaRotorDiskSource::setFaceArea(vector& axis, const bool correct)
{
    area_ = 0.0;
    areaSf_ = vector::zero;

    static const scalar tol = 0.999; //0.8  CAE/RGK - tune this

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = Zero;

    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    UIndirectList<label>(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nFaces() - nInternalFaces, -1);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // Add internal field contributions
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[facei];
                areaSf_[own] += Sf[faceI];
                n += Sf[facei];
                rotorDiskV_ += mesh_.V()[own]; // CAE/RGK 2017/09/13
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[facei];
                areaSf_[nbr] -= Sf[faceI];
                n -= Sf[facei];
            }
        }
    }


    // Add boundary contributions
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const label nbr = nbrFaceCellAddr[facei - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }

    if (correct)
    {
        reduce(n, sumOp<vector>());
        axis = n/mag(n);
    }

    if (debug)
    {
        volScalarField area
        (
            IOobject
            (
                name_ + ":area",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimArea, 0)
        );
        UIndirectList<scalar>(area.primitiveField(), cells_) = area_;

        Info<< type() << ": " << name_ << " writing field " << area.name()
            << endl;

        area.write();
    }

    // CAE/RGK 2017/09/13
    reduce(rotorDiskV_, sumOp<scalar>());

    Info<< "Volume of rotor disk cells = " << rotorDiskV_ << endl;



}


void Foam::fv::panaRotorDiskSource::createCoordinateSystem()
{
    // Construct the local rotor co-prdinate system
    vector origin(Zero);
    vector axis(Zero);
    vector refDir(Zero);

    geometryModeType gm =
        geometryModeTypeNames_.read(coeffs_.lookup("geometryMode"));

    switch (gm)
    {
        case gmAuto:
        {
            // Determine rotation origin (cell volume weighted)
            scalar sumV = 0.0;
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                sumV += V[celli];
                origin += V[celli]*C[celli];
            }
            reduce(origin, sumOp<vector>());
            reduce(sumV, sumOp<scalar>());
            origin /= sumV;

            // Determine first radial vector
            vector dx1(Zero);
            scalar magR = -GREAT;
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                vector test = C[celli] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }
            }
            reduce(dx1, maxMagSqrOp<vector>());
            magR = mag(dx1);

            // Determine second radial vector and cross to determine axis
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                vector dx2 = C[celli] - origin;
                if (mag(dx2) > 0.5*magR)
                {
                    axis = dx1 ^ dx2;
                    if (mag(axis) > SMALL)
                    {
                        break;
                    }
                }
            }
            reduce(axis, maxMagSqrOp<vector>());
            axis /= mag(axis);

            // Correct the axis direction using a point above the rotor
            {
                vector pointAbove(coeffs_.lookup("pointAbove"));
                vector dir = pointAbove - origin;
                dir /= mag(dir);
                if ((dir & axis) < 0)
                {
                    axis *= -1.0;
                }
            }

            coeffs_.lookup("refDirection") >> refDir;

            // set the face areas and apply correction to calculated axis
            // e.g. if cellZone is more than a single layer in thickness
            setFaceArea(axis, true);

            break;
        }
        case gmSpecified:
        {
            coeffs_.lookup("origin") >> origin;
            coeffs_.lookup("axis") >> axis;
            axis /=  mag(axis);

            //coeffs_.lookup("refDirection") >> refDir;
            // Determine a reference direction which lies normal
            // to both the axis and any radial vector.

            // Determine face areas in faceZone
            setFaceArea(axis, false);

            // Determine any radial vector
             const vectorField& C = mesh_.C();
            vector dx1(vector::zero);
            rotorDiskA_ = gSum( area_ ); // CAE/RGK 2017/09/13
            // Employ a finite size mag(test_radial_vector)
            scalar magR = 0.1 * sqrt( rotorDiskA_ / M_PI);

            // Extract test radial vector
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                vector test = C[cellI] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }

            }
            reduce(dx1, maxMagSqrOp<vector>());
            //magR = mag(dx1);

            refDir = axis ^ dx1;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown geometryMode " << geometryModeTypeNames_[gm]
                << ". Available geometry modes include "
                << geometryModeTypeNames_ << exit(FatalError);
        }
    }


    // CAE/RGK 2017/09/13 - use the sign of omega to determine the sense of the rotation axis
    // so that we only consider anti-clockwise rotations.
    scalar rpm(readScalar(coeffs_.lookup("rpm")));
    axis *=  Foam::sign(rpm);

    coordSys_ = cylindricalCS("rotorCoordSys", origin, axis, refDir, false);
    const scalar diameter = Foam::sqrt(4.0*rotorDiskA_/mathematical::pi);
    Info<< "    Rotor gometry:" << nl
        << "    - disk diameter = " << diameter << nl
        << "    - disk area     = " << rotorDiskA_ << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.R().e1() << nl
        << "    - psi-axis      = " << coordSys_.R().e2() << nl
        << "    - refDir    = " << refDir << nl
        << "    - z-axis        = " << coordSys_.R().e3() << endl;

}


void Foam::fv::panaRotorDiskSource::constructGeometry()
{
    const vectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {

            const label cellI = cells_[i];

            // position in (planar) rotor co-ordinate system
            x_[i] = coordSys_.localPosition(C[cellI]);

            rMax_ = max(rMax_, x_[i].x());

            // swept angle relative to rDir axis [radians] in range 0 -> 2*pi
            scalar psi = x_[i].y();

            // blade flap angle [radians]
            scalar beta =
                flap_.beta0 - flap_.beta1c*cos(psi) - flap_.beta2s*sin(psi);

            // determine rotation tensor to convert from planar system into the
            // rotor cone system
            scalar c = cos(beta);
            scalar s = sin(beta);
            R_[i] = tensor(c, 0, -s, 0, 1, 0, s, 0, c);
            invR_[i] = R_[i].T();

            // CAE/RGK - 2017/08 - addition
            // Apply cylindrical transformation

            c = cos(psi);
            s = sin(psi);
            Rcyl_[i] = tensor(c, s, 0, -s, c, 0, 0, 0, 1);
            invRcyl_[i] = Rcyl_[i].T();

            //Pout<< "  R_[i] = " <<  R_[i] << endl;
            //Pout<< "  invR_[i] = " <<  invR_[i] << endl;

        }// if (area_[i] > ROOTVSMALL)

    }//forAll(cells_, i)

}


Foam::tmp<Foam::vectorField> Foam::fv::panaRotorDiskSource::inflowVelocity
(
    const volVectorField& U
) const
{
    switch (inletFlow_)
    {
        case ifFixed:
        case ifSurfaceNormal:
        {
            return tmp<vectorField>
            (
                new vectorField(mesh_.nCells(), inletVelocity_)
            );

            break;
        }
        case ifLocal:
        {
            return U.internalField();

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown inlet flow specification" << abort(FatalError);
        }
    }

    return tmp<vectorField>(new vectorField(mesh_.nCells(), vector::zero));
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::panaRotorDiskSource::panaRotorDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr

)
:
    cellSetOption(name, modelType, dict, obr),
    rhoName_("none"),
    rhoRef_(1.0),
    omega_(0.0),
    nBlades_(0),
    inletFlow_(ifLocal),
    flowRate0_(0.0),
    flowRate_(0.0),
    fanAxis_(vector::zero),
    inletVelocity_(vector::zero),
    feedbackTime_(0),
    tipEffect_(1.0),
    flap_(),
    x_(cells_.size(), vector::zero),
    R_(cells_.size(), I),
    Rcyl_(cells_.size(), I),
    invR_(cells_.size(), I),
    invRcyl_(cells_.size(), I),
    area_(cells_.size(), 0.0),
    areaSf_(cells_.size(), vector::zero),
    rotorDiskA_(0.0),
    coordSys_(false),
    rMax_(0.0),
    rotorDiskV_(0.0),
    trim_(trimModel::New(*this, coeffs_)),
    blade_(coeffs_.subDict("blade")),
    profiles_(coeffs_.subDict("profiles"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::panaRotorDiskSource::~panaRotorDiskSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::panaRotorDiskSource::calculate
(
    const vectorField& U,
    const scalarField& thetag,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh_.V();
    const bool compressible = this->compressible();
    tmp<volScalarField> trho(rho());

    // logging info
    scalar dragEff = 0.0;
    scalar liftEff = 0.0;
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;

    forAll(cells_, i)
    {

        const label cellI = cells_[i];

        if (area_[i] > ROOTVSMALL)
        {
            //const label cellI = cells_[i];

            const scalar radius = x_[i].x();

            // velocity in local cylindrical reference frame
            //vector Uc = coordSys_.localVector(U[cellI]);
            vector Uc = coordSys_.R().R().T() & U[cellI];

            // transform from rotor cylindrical into local coning system
            Uc = R_[i] & Uc;
            Uc = Rcyl_[i] & Uc;

            // set radial component of velocity to zero
            Uc.x() = 0.0;

            // set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();

            // determine blade data for this radius
            // i2 = index of upper radius bound data point in blade list
            scalar twist = 0.0;
            scalar chord = 0.0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0.0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);

            // flip geometric angle if blade is spinning in reverse (clockwise)
            scalar alphaGeom = thetag[i] + twist;
            if (omega_ < 0)
            {
                alphaGeom = mathematical::pi - alphaGeom;
            }

            // effective angle of attack
            scalar alphaEff = alphaGeom - atan2(-Uc.z(), Uc.y());

            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // determine profile data for this radius and angle of attack
            const label profile1 = blade_.profileID()[i1];
            const label profile2 = blade_.profileID()[i2];

            scalar Cd1 = 0.0;
            scalar Cl1 = 0.0;
            profiles_[profile1].Cdl(alphaEff, Cd1, Cl1);

            scalar Cd2 = 0.0;
            scalar Cl2 = 0.0;
            profiles_[profile2].Cdl(alphaEff, Cd2, Cl2);

            scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // apply tip effect for blade lift
            scalar tipFactor = neg(radius/rMax_ - tipEffect_);

            // calculate forces perpendicular to blade
            //scalar pDyn = 0.5*magSqr(Uc);
            // CAE/RGK  - dynamic pressure should be obtained from rotational motion.
            scalar pDyn = 0.5*magSqr(Uc.y());

            if (compressible)
            {
                pDyn *= trho()[cellI];
            }

            scalar f = pDyn*chord*nBlades_*area_[i]/radius/mathematical::twoPi;
            // CAE/RGK - specify Cl = 0 in fvOptions. Normal flow is represented by a specified
            // flow rate.
            vector localForce = vector(0.0, -f*Cd, tipFactor*f*Cl);

            // accumulate forces
            dragEff += rhoRef_*localForce.y();
            liftEff += rhoRef_*localForce.z();

            // convert force from local coning system into rotor cylindrical
            localForce = invRcyl_[i] & localForce;
            localForce = invR_[i] & localForce;

            // convert force to global cartesian co-ordinate system
            // force[cellI] = coordSys_.globalVector(localForce); CAE/RGK - error, returns 0. Corrected in OF4.x?
            force[cellI] =   (coordSys_.R().R() & localForce);

            if (compressible)
            {
                force[cellI]  *= trho()[cellI];
            }

            if (divideVolume)
            {
                force[cellI] /= V[cellI];
            }

        }//   if (area_[i] > ROOTVSMALL)

        // CAE/RGK:  Subtract constant momentum flux term per unit vol. (this is included as a
        // negative source term, see below at line 730).
        // V_ : total volume of cells in cellZone
        // force[cellI]  -=   inletVelocity_ * flowRate_ / V_;
        force[cellI]  -=   fanAxis_ * flowRate_  * flowRate_  / ( rotorDiskA_  * V_);

    }// forAll(cells_, i)

    if (output)
    {
        reduce(AOAmin, minOp<scalar>());
        reduce(AOAmax, maxOp<scalar>());
        reduce(dragEff, sumOp<scalar>());
        reduce(liftEff, sumOp<scalar>());

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << "    Effective drag = " << dragEff << nl
            << "    Effective lift = " << liftEff << endl;

    }

}


void Foam::fv::panaRotorDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    dimensionSet dims = dimless;
    if (eqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoName") >> rhoName_;
        dims.reset(dimForce/dimVolume);
    }
    else
    {
        coeffs_.lookup("rhoRef") >> rhoRef_;
        dims.reset(dimForce/dimVolume/dimDensity);
    }

    volVectorField force
    (
        IOobject
        (
            name_ + ":panaRotorForce",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dims, vector::zero)
    );

    const volVectorField& U = eqn.psi();

    const vectorField Uin(inflowVelocity(U));

    trim_->correct(Uin, force);

    calculate(Uin, trim_->thetag(), force);

    // add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
}


void Foam::fv::panaRotorDiskSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::panaRotorDiskSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {

        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // read co-ordinate system/geometry invariant properties
        scalar rpm(readScalar(coeffs_.lookup("rpm")));
        omega_ = rpm/60.0*mathematical::twoPi;

        Info<< "- rpm = " << rpm << endl;

        coeffs_.lookup("nBlades") >> nBlades_;

        Info<< "Number of blades = " << nBlades_ << endl;

        inletFlow_ = inletFlowTypeNames_.read(coeffs_.lookup("inletFlowType"));

        feedbackTime_ = coeffs_.lookupOrDefault<scalar>("feedbackTime", 1000);


        coeffs_.lookup("tipEffect") >> tipEffect_;

        const dictionary& flapCoeffs(coeffs_.subDict("flapCoeffs"));
        flapCoeffs.lookup("beta0") >> flap_.beta0;
        flapCoeffs.lookup("beta1c") >> flap_.beta1c;
        flapCoeffs.lookup("beta2s") >> flap_.beta2s;
        flap_.beta0 = degToRad(flap_.beta0);
        flap_.beta1c = degToRad(flap_.beta1c);
        flap_.beta2s = degToRad(flap_.beta2s);

        // create co-ordinate system
        createCoordinateSystem();

        // read co-odinate system dependent properties
        checkData();

        constructGeometry();

        trim_->read(coeffs_);

        if (debug)
        {
            writeField("thetag", trim_->thetag()(), true);
            writeField("faceArea", area_, true);
        }

        return true;
    }
    else
    {
        return false;
    }
}


//- Scalar CAE/RGK 2017/09/28
void Foam::fv::panaRotorDiskSource::correct(volVectorField& fld)
{

    // CAE/RGK 2017/10/02
    // Confirm flux emanating from fan
    scalar  flowRateLocal = 0.0;

    forAll(cells_, i)
    {

        if (area_[i] > ROOTVSMALL)
        {
            const label cellI = cells_[i];
            flowRateLocal += (  fld.internalField()[cellI]  & areaSf_[i] );
        }

    }

    reduce(flowRateLocal, sumOp<scalar>());

    Info<< "panaRotorDisk::correct():: Flux computed at fan  = "<< flowRateLocal
     << "   ( user-specified = "<<flowRate0_ <<" )"<< endl;

    // CAE/RGK 2017/10/02
    // Apply feedback

    if (mesh_.time().timeOutputValue() > feedbackTime_)
    {
            flowRate_  += ( flowRate0_ - flowRateLocal);
    }

}

// ************************************************************************* //
