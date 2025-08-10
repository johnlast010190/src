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
    (c) 2023 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "radialThrustActuationDiskSource.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvm/fvm.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::radialThrustActuationDiskSource::
addRadialActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    sampleDisk();
    wakeFraction();
    controlDisk();

    if (shipResistanceThrust_)
    {
        functionObjects::forces shipForce(word("forces"), mesh(), coeffs_.subDict("forces"));
        shipForce.calcForcesMoment();
        const_cast<radialThrustActuationDiskSource*>(this)
            ->updateThrust(shipForce.forceEff().x()-SFC_);
    }

    scalarField Tr(cells.size());
    const vector uniDiskDir = diskDir_/mag(diskDir_);

    const Field<vector> zoneCellCentres(mesh().cellCentres(), cells);
    const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), cells);

    scalar a = radialCoeffs_[0];
    scalar m = radialCoeffs_[1];
    scalar n = radialCoeffs_[2];
    vector outputThrust = vector(0, 0, 0);
    scalar torque = 0;
    vector testThrust    = vector(0, 0, 0);
    scalar A  = 1.0;
    scalar f1 = 0;

    // determine the new position of the center of the zone
    vector newPropPosition = vector(0,0,0);
    forAll(cells, i)
    {
        newPropPosition +=  mesh().cellCentres()[cells[i]]*Vcells[cells[i]]/V_;
    }

    reduce(newPropPosition, sumOp<vector>());
    // Calculate scale factor A
    forAll(cells, i)
    {
        vector radRel = propOrientation_.invTransform
        (
            mesh().cellCentres()[cells[i]] - newPropPosition
        );
        scalar radius =
            sqrt(sqr(radRel.component(1)) + sqr(radRel.component(2)));

        if (radius < R_ && radius > rh_)
        {
            scalar radCarat    = (radius - rh_)/(R_ - rh_ +VSMALL);
            f1 = pow(mag(radCarat), m)* pow( mag(a - radCarat)/(a + VSMALL), n);
        }
        testThrust += (Vcells[cells[i]]/V_)*f1*uniDiskDir;
        f1 = 0;
    }
    reduce(testThrust, sumOp<vector>());
    A  = 1.0/(mag(testThrust) + VSMALL);

    f1 = 0;
    forAll(cells, i)
    {
        vector radRel = propOrientation_.invTransform
        (
            mesh().cellCentres()[cells[i]] - newPropPosition
        );
        scalar radius =
            sqrt(sqr(radRel.component(1)) + sqr(radRel.component(2)));
        scalar theta = atan2( radRel.component(2), radRel.component(1) );

        if (radius < R_ && radius > rh_)
        {
            scalar radCarat = (radius - rh_)/(R_ - rh_ +VSMALL);
            f1 = T()*A*pow(mag(radCarat), m)
                * pow( mag(a - radCarat)/(a + VSMALL), n);
        }

        vector fz_ = vector(0, 0, 1)*PoverD_*0.3183*(R_/(radius + VSMALL))*f1;
        vector fy_ = vector(0, 1, 0)*PoverD_*-0.3183*(R_/(radius + VSMALL))*f1;

        Usource[cells[i]] += (Vcells[cells[i]]/V_)*
            propOrientation_.transform( f1*uniDiskDir + fz_*cos(theta) + fy_*sin(theta) );

        outputThrust += (Vcells[cells[i]]/V_)*
            propOrientation_.transform( f1*uniDiskDir + fz_*cos(theta) + fy_*sin(theta) );

        torque += (Vcells[cells[i]]/V_)
            *f1*PoverD_*0.3183*(R_/(radius + VSMALL))*radius;

        f1 = 0;
    }
    reduce(torque, sumOp<scalar>());
    reduce(outputThrust, sumOp<vector>());

    Info<< "Input Thrust " << T() << nl
        << "Force factor " << A  << nl
        << "Output Thrust " << outputThrust.component(0) << nl
        << "Torque " << torque << nl
        << "Propeller Centre" << propPosition_ << nl
        << "NewPropPosition " << newPropPosition << endl;

}


// ************************************************************************* //
