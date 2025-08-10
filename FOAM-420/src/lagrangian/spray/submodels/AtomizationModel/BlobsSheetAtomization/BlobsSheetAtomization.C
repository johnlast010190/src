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
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/AtomizationModel/BlobsSheetAtomization/BlobsSheetAtomization.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BlobsSheetAtomization<CloudType>::BlobsSheetAtomization
(
    const dictionary& dict,
    CloudType& owner
)
:
    AtomizationModel<CloudType>(dict, owner, typeName),
    B_(readScalar(this->coeffDict().lookup("B"))),
    angle_(readScalar(this->coeffDict().lookup("angle")))
{}


template<class CloudType>
Foam::BlobsSheetAtomization<CloudType>::BlobsSheetAtomization
(
    const BlobsSheetAtomization<CloudType>& am
)
:
    AtomizationModel<CloudType>(am),
    B_(am.B_),
    angle_(am.angle_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BlobsSheetAtomization<CloudType>::~BlobsSheetAtomization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::BlobsSheetAtomization<CloudType>::initLiquidCore() const
{
    return 1.0;
}


template<class CloudType>
bool Foam::BlobsSheetAtomization<CloudType>::calcChi() const
{
    return false;
}


template<class CloudType>
void Foam::BlobsSheetAtomization<CloudType>::update
(
    const scalar dt,
    scalar& d,
    scalar& liquidCore,
    scalar& tc,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const scalar volFlowRate,
    const scalar rhoAv,
    const scalar Urel,
    const vector& pos,
    const vector& injectionPos,
    const scalar pAmbient,
    const scalar chi,
    Random& rndGen
) const
{
    scalar lBU =
        B_
      * sqrt
        (
            rho*sigma*d*cos(angle_*constant::mathematical::pi/360.0)
          / sqr(rhoAv*Urel)
        );

    scalar pWalk = mag(pos - injectionPos);

    if (pWalk > lBU)
    {
        liquidCore = 0.0;
    }
}


// ************************************************************************* //
