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
    (c) 2010-2022 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "thrusterActuationDiskSource.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "finiteVolume/fvm/fvm.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::thrusterActuationDiskSource::
addThrusterActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    scalarField dP(cells.size(), 0.0);

    vector aveVel = vector::zero;

    switch (fanModel_)
    {
        case fsAverageVelocity:
        {
            // Enter the dP vs U curve with the average velocity
            forAll(cells, i)
            {
                const label celli = cells[i];
                aveVel += Vcells[celli]*U[celli];
            }
            reduce(aveVel, sumOp<vector>());
            aveVel /= V_;

            const scalar Un =
                min
                (
                    max
                    (
                        (aveVel & coorFramePtr_->axis()),
                        lowerBound_
                    ),
                    upperBound_
                );
            dP = f_.value(Un);

            // Reporting working point
            Info<< indent << "Report thruster source: "
                << this->name() << nl
                << indent << "dP [Pa] = " << dP << nl
                << indent << "Uave [m/s] = " << Un << nl
                << indent << nl;
            break;
        }

        case fsLocalVelocity:
        {
            // Enter the dP vs U curve with the local velocity
            forAll(cells, i)
            {
                const label celli = cells[i];
                const scalar Un =
                    min
                    (
                        max
                        (
                            (U[celli] & coorFramePtr_->axis()),
                            lowerBound_
                        ),
                        upperBound_
                    );
                dP[i] = f_.value(Un);
            }

            // Reporting working point
            const scalar flowRate =
                0.5*this->zoneFlux(mesh_, zoneBoundaryFaces_);
            const scalar Uave = flowRate/diskArea_;

            Info<< indent << "Report thruster source: "
                << this->name() << nl
                << indent << "dP [Pa] = " << f_.value(Uave) << nl
                << indent << "Uave [m/s] = " << Uave << nl
                << indent << nl;
            break;
        }

        case fsFlowRate:
        {
            // enter the dP vs Q with the flowrate across the cell zone
            const scalar flowRate =
                min
                (
                    max
                    (
                        (0.5*this->zoneFlux(mesh_, zoneBoundaryFaces_)),
                        lowerBound_
                    ),
                    upperBound_
                );
            dP = f_.value(flowRate);

            Info<< indent << "Report thruster source: "
                << this->name() << nl
                << indent << "dP [Pa] = " << dP << nl;

            if (!compressible_)
            {
                Info<< indent << "Flow Rate [m3/s] = " << flowRate << nl;
            }
            else
            {
                Info<< indent << "Flow Rate [kg/s] = " <<flowRate << nl;
            }
            Info<< indent << nl;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown fan source type. Valid types are: "
                << fanModelTypeNames_ << nl << exit(FatalError);
        }
    }

    forAll(cells, i)
    {
        const label celli = cells[i];
        Usource[celli] +=
            Vcells[celli]
           *(
                (-load_[i]*dP[i]/rhoRef_*coorFramePtr_->axis())
              + rho[celli]*jumpV_[i]*(U[celli] & coorFramePtr_->axis())
            );
    }

}


// ************************************************************************* //
