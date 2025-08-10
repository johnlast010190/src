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
    (c) 2018 Esi Ltd.
    (c) 2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::powerLawAnisotropic::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{

    const scalar  Bm1over2 = (B_ - 1.0)/2.0;

    vectorField Urel = U;

    if (mesh_.moving() && porousRelVelocity_)
    {
        //mesh is moving
        Urel -= fvc::reconstruct
        (
            (fvc::meshPhi(mesh_.lookupObject<volVectorField>("U")))
        )->internalField();
    }

    if (coorFramePtr_)
    {
        Urel -= coorFramePtr_->frameVelocity(mesh_.C().internalField(), true);
    }

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];

            const tensor Cd = rho[celli]*C_*pow(magSqr(Urel[celli]), Bm1over2);

            const scalar isoCd = tr(Cd);

            Udiag[celli] += V[celli]*isoCd;
            Usource[celli] -= V[celli]*((Cd & Urel[celli])  - ((I*isoCd) & U[celli]));
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLawAnisotropic::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    if (validParentFrame())
    {
        FatalErrorInFunction
            << "Tensorial anisotropic powerLaw not supported in moving "
            << "reference frame." << exit(FatalError);
    }

    const scalar  Bm1over2 = (B_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            AU[celli] += rho[celli]*C_*pow(magSqr(U[celli]), Bm1over2);
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLawAnisotropic::adjointApply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& Uprimal,
    const vectorField& U
) const
{

    vectorField UpRel = Uprimal;
    if (coorFramePtr_)
    {
        UpRel -= coorFramePtr_->frameVelocity(mesh_.C().internalField(), true);
    }
    const scalar  Bm1over2 = (B_ - 1.0)/2.0;
    const tensor& Ct = C_.T();

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
        forAll(cells, i)
        {
            const vector& Upr(UpRel[cells[i]]);

            const tensor dragCoeff
                = rho[cells[i]]*Ct*pow(magSqr(Upr), Bm1over2);

            const scalar isoDragCoeff = tr(dragCoeff);

            Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
            Usource[cells[i]] -=
                (V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]));

        }
    }
}

// ************************************************************************* //
