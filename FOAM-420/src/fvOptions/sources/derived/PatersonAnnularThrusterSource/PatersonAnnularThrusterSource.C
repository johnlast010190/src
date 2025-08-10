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
    (c) 2016-2022 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "PatersonAnnularThrusterSource.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "sets/cellSources/cylinderAnnulusToCell/cylinderAnnulusToCell.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(PatersonAnnularThrusterSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        PatersonAnnularThrusterSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::PatersonAnnularThrusterSource::PatersonAnnularThrusterSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    shipForce_(word("shipForce"), mesh_, coeffs_.subDict("shipForce")),
    Uship_(coeffs_.lookup("Uship")),
    rhoRef_(readScalar(coeffs_.subDict("shipForce").lookup("rhoInf"))),
    R_(readScalar(coeffs_.lookup("propTipRadius"))),
    rh_(readScalar(coeffs_.lookup("propHubRadius"))),
    propPosition_(coeffs_.lookup("propPosition")),
    diskDir_(coeffs_.lookup("propDirection")),
    tProp_(readScalar(coeffs_.lookup("propThickness"))),
    qPropDir_(coeffs_.lookupOrDefault<scalar>("torqueDirection", 1.0)),
    diskFraction_(readScalar(coeffs_.lookup("diskFraction"))),
    J_(readScalar(coeffs_.lookup("J"))),
    a_(coeffs_.lookup("a")),
    b_(coeffs_.lookup("b")),
    normaliseForces_(coeffs_.lookupOrDefault<Switch>("normaliseForces", true)),
    wakeFraction_(coeffs_.lookupOrDefault<Switch>("wakeFraction", false))
{
    Info<< "    - creating equilibrium actuation disk: " << name_ << endl;

    fieldNames_ = List<word>(1, "U");
    applied_.setSize(fieldNames_.size(), false);

    diskDir_ /= mag(diskDir_);
    checkData();

    //select the cells in the annular disk
    setCellSet();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::fv::PatersonAnnularThrusterSource::calcDrag()
{
    shipForce_.calcForcesMoment();

    return
    (
        shipForce_.forceEff()
        /diskFraction_
    );
}

Foam::vector Foam::fv::PatersonAnnularThrusterSource::calcTorque
(
    const vector& thrust
)
{
    scalar magThrust(mag(thrust));
    scalar magUship(mag(Uship_));

    if (magThrust < SMALL || magUship < SMALL)
    {
        WarningInFunction
            << "Small ship velocity or thrust detected (Uship = "
            << magUship << ", thrust = " << magThrust
            << "). Setting torque to zero." << endl;

        return vector::zero;
    }


    // iteration variables
    scalar tol(1e-6);
    label maxIter(100);
    scalar J0 = J_;
    label iter = 0;

    // set initial guess for rps (n) and Kt
    scalar Dp(2*R_);
    scalar n = magUship/(J0*Dp);
    scalar thrustDivRhoU2Dp2(magThrust/(rhoRef_*pow(magUship,2)*pow(Dp,2)));

    // Iteration loop
    do
    {
        iter++;

        // Evaluate the function and its derivative
        // f is a 5th-order polynomial curve-fit for KT(J)
        scalar f = a_.value(J_) - pow(J_,2)*thrustDivRhoU2Dp2;
        scalar dfdx = a_.derivative(J_) - 2.*J0*thrustDivRhoU2Dp2;

        // Newton-Raphson step
        J0 = J_;

        if (iter > maxIter/2.0) //relax update
        {
            J_ -= 0.5*f/dfdx;
        }
        else
        {
            J_ -= f/dfdx;
        }

    } while (mag(J_ - J0)/mag(J_) > tol && iter < maxIter);

    n = magUship/(J_*Dp);
    //scalar Kt(magThrust / (rhoRef_ * pow(n,2) * pow(Dp,4)));
    scalar Kq = 0.1 * b_.value(J_); //whatr is the "0.1" for?
    //scalar Kq = b_.value(J_);

    //removed "pos0" function to allow negative rotation around thrust axis
    //useful for dual propeller scenarious
    scalar torque(pos0(Kq)*Kq * rhoRef_ * pow(n,2) * pow(Dp,5));

    return qPropDir_ * torque * thrust / magThrust;
}

Foam::scalar Foam::fv::PatersonAnnularThrusterSource::annularVolume() const
{
    return
    (
        constant::mathematical::pi * (R_*R_ - rh_*rh_)*tProp_
    );
}

void Foam::fv::PatersonAnnularThrusterSource::checkData() const
{
    if (annularVolume() <= VSMALL)
    {
        FatalErrorInFunction
           << "disk volume is approximately zero"
           << exit(FatalIOError);
    }
}
void Foam::fv::PatersonAnnularThrusterSource::setCellSet()
{
    cylinderAnnulusToCell propShape
    (
        mesh_,
        propPosition_,
        propPosition_ + diskDir_*tProp_,
        R_,
        rh_
    );

    cellSet propellerCellSet
    (
        mesh_,
        "propellerCellSet",
        mesh_.nCells()/10+1  // Reasonable size estimate.
    );

    propShape.applyToSet(topoSetSource::NEW, propellerCellSet);

    cells_ = propellerCellSet.toc();

    V_ = 0.0;
    forAll(cells_, i)
    {
        V_ += mesh_.V()[cells_[i]];
    }
    reduce(V_, sumOp<scalar>());

    if (!mesh_.changing())
    {
        Info<< "- selected " << returnReduce(cells_.size(), sumOp<label>())
            << " cell(s) with volume " << V_ << nl << endl;
    }
}

void Foam::fv::PatersonAnnularThrusterSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        addBodyForce
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U
        );
    }
}

void Foam::fv::PatersonAnnularThrusterSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        addBodyForce
        (
            Usource,
            cells_,
            cellsV,
            rho,
            U
        );
    }
}


bool Foam::fv::PatersonAnnularThrusterSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        shipForce_.read(coeffs_.subDict("shipForce"));
        coeffs_.lookup("Uship")>> Uship_;
        coeffs_.subDict("shipForce").readIfPresent("rhoInf", rhoRef_);
        coeffs_.readIfPresent("R", R_);
        coeffs_.readIfPresent("rh", rh_);
        coeffs_.readIfPresent("propPosition", propPosition_);
        coeffs_.readIfPresent("propDirection", diskDir_);
        coeffs_.readIfPresent("propThickness", tProp_);
        qPropDir_ = coeffs_.lookupOrDefault<scalar>("torqueDirection", 1.0);
        diskFraction_ = coeffs_.lookupOrDefault<scalar>("diskFraction", 1.0);
        coeffs_.readIfPresent("J", J_);
        normaliseForces_
            = coeffs_.lookupOrDefault<Switch>("normaliseForces", true);
        wakeFraction_
            = coeffs_.lookupOrDefault<Switch>("wakeFraction", false);

        a_ = Polynomial<6>(coeffs_.lookup("a"));
        b_ = Polynomial<6>(coeffs_.lookup("b"));

        diskDir_ /= mag(diskDir_);
        checkData();

        //select the cells in the annular disk
        setCellSet();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
