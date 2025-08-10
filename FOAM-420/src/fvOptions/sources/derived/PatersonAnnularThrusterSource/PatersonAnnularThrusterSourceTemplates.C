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
    (c) Esi Ltd

\*----------------------------------------------------------------------------*/

#include "PatersonAnnularThrusterSource.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "finiteVolume/fvm/fvm.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::PatersonAnnularThrusterSource::addBodyForce
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
)
{
    //use function object to calculate drag
    vector drag(calcDrag());

    //drag in ship propagation direction (needed for half-model or unconverged)
    vector shipDir(Uship_/mag(Uship_));
    drag = max(0, (drag & shipDir)) * shipDir;

    //project drag onto actuator disk normal
    vector thrust = drag / (shipDir & diskDir_);

    vector torque(calcTorque(thrust));

    scalar integratedThrust(0.0);
    scalar integratedTorque(0.0);

    if (V_ > 0)
    {
        List<vector> thrustSource(cells.size(), vector::zero);
        List<vector> torqueSource(cells.size(), vector::zero);
        scalar rVfrac = diskFraction_/V_;

        forAll(cells, i)
        {
            label cellI = cells[i];

            // dimensionless radius
            vector cellRadius(mesh_.C()[cellI] - propPosition_);
            cellRadius -= diskDir_ * (cellRadius & diskDir_);

            scalar magr(mag(cellRadius));
            scalar rstar
            (
                max(0, (magr - rh_))/(R_ - rh_)
            );

            scalar thrustDist = 3.5183 * rstar * ::sqrt(1.0 - rstar);

            scalar torqueDist = 2.3757 * rstar * ::sqrt(1.0-rstar)
                       /((1.0-rh_/R_)*rstar + rh_/R_);

            //axial cell force
            vector dThrust
                = thrust* thrustDist * rVfrac * mesh_.V()[cellI];

            // tangential cell force magnitude
            vector dTorque
                = torque * torqueDist* rVfrac * mesh_.V()[cellI] / magr;

            thrustSource[i] = dThrust;
            //correct torque force direction
            torqueSource[i] = (dTorque ^ (cellRadius/magr));

            // for diagnostic and scaling, compute total integrated
            //thrust & torque
            integratedThrust += mag(dThrust);
            integratedTorque += mag(dTorque)*magr; //force times distance
        }

        reduce(integratedThrust,sumOp<scalar>());
        reduce(integratedTorque,sumOp<scalar>());

        if (normaliseForces_)
        {
            scalar tScale(mag(thrust)/max(SMALL, integratedThrust));
            scalar qScale(mag(torque)/max(SMALL, integratedTorque));

            thrustSource = tScale * thrustSource;
            torqueSource = qScale * torqueSource;
        }

        if (isA<geometricOneField>(rho))
        {
            forAll(cells, i)
            {
                label cellI = cells[i];

                Usource[cellI]
                    -= (thrustSource[i] + torqueSource[i])/rhoRef_;
            }
        }
        else
        {
            forAll(cells, i)
            {
                label cellI = cells[i];
                Usource[cellI]
                    -= thrustSource[i] + torqueSource[i];
            }

        }

        if (wakeFraction_)
        {
            // compute Va (= zonal_Umean) and wake fraction
            vector Va(vector::zero);
            scalar Umin(VGREAT);
            scalar Umax(VSMALL);
            scalar w(0);

            forAll(cells, i)
            {
                label cellI = cells[i];

                // compute Va
                Va += U[cellI] * mesh_.V()[cellI];
                if (Umin > mag(U[cellI])) Umin = mag(U[cellI]);
                if (Umax < mag(U[cellI])) Umax = mag(U[cellI]);
            }

            reduce(Va,sumOp<vector>());
            scalar Vglob = V_;
            reduce(Vglob,sumOp<scalar>());
            reduce(Umin,minOp<scalar>());
            reduce(Umax,maxOp<scalar>());

            Va /= Vglob;
            w = ((Uship_ - Va) & Uship_)/sqr(mag(Uship_));

            Info<< "advance speed  = " << Va << endl;
            Info<< "wake fraction  = " << w << endl;
            Info<< "mag(Va)        = " << mag(Va) << endl;
            Info<< "min(V), max(V) = " << Umin << " , " << Umax << endl;
        }

        if (debug)
        {

            Info<< "thrust                        = " << thrust << endl;
            Info<< "torque                        = " << torque << endl;
            Info<< "integratedThrust              = "
                 << integratedThrust << endl;
            Info<< "integratedTorque              = "
                 << integratedTorque << endl;
            Info<< "integratedThrust/diskFraction = "
                 << integratedThrust/diskFraction_ << endl;
            Info<< "integratedTorque/diskFraction = "
                 << integratedTorque/diskFraction_ << endl;

            if (normaliseForces_)
            {
                Info<< "Automatic force normalisation enabled." << endl;
                Info<< "normalised integratedThrust force = "
                     << gSum(thrustSource) << endl;

                vector torqueMoment(vector::zero);
                forAll(cells, i)
                {
                    label cellI = cells[i];
                    vector cellRadius(mesh_.C()[cellI] - propPosition_);
                    cellRadius -= diskDir_ * (cellRadius & diskDir_);
                    torqueMoment += (cellRadius ^ torqueSource[i]);
                }
                reduce(torqueMoment,sumOp<vector>());
                Info<< "normalised integratedTorque moment) = "
                     << torqueMoment << endl;
            }
        }
    }


}


// ************************************************************************* //
