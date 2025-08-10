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
    (c) 2010-2016 Esi Ltd.
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::powerLaw::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar t = db().time().timeOutputValue();
    const scalar C0 = C0_->value(t);
    const scalar C1m1b2 = (C1_->value(t) - 1.0)/2.0;

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

            scalar Pcoeff = V[celli]*rho[celli]*C0
                *pow(magSqr(Urel[celli]), C1m1b2);
            Udiag[celli] += Pcoeff;
            Usource[celli] += Pcoeff * (U[celli] - Urel[celli]);
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLaw::apply
(
    tensorField& AU,
    vectorField& source,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    vectorField Urel = U;

    if (mesh_.moving() && porousRelVelocity_)
    {
        //mesh is moving
        Urel -= fvc::reconstruct
        (
            (fvc::meshPhi(mesh_.lookupObject<volVectorField>("U")))
        )->internalField();
    }

    tmp<vectorField> tframeU;

    if (coorFramePtr_)
    {
        tframeU =
            coorFramePtr_->frameVelocity(mesh_.C().internalField(), true);
        Urel -= tframeU();
    }

    const scalar t = db().time().timeOutputValue();
    const scalar C0 = C0_->value(t);
    const scalar C1m1b2 = (C1_->value(t) - 1.0)/2.0;
    const scalarField& V = mesh_.V();

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            scalar Pcoeff = V[celli]*rho[celli]*C0*
                pow(magSqr(Urel[celli]), C1m1b2);

            tensor tPcoeff = I*Pcoeff;

            AU[celli] += tPcoeff;

            if (coorFramePtr_)
            {
                source[celli] += tPcoeff&tframeU()[celli];
            }
        }
    }
}

template<class RhoFieldType>
void Foam::porosityModels::powerLaw::adjointApply
(
    scalarField& Udiag,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& Uprimal
) const
{
    vectorField UpRel = Uprimal;
    if (coorFramePtr_)
    {
        UpRel -= coorFramePtr_->frameVelocity(mesh_.C().internalField(), true);
    }
    const scalar t = db().time().timeOutputValue();
    const scalar C0 = C0_->value(t);
    const scalar C1m1b2 = (C1_->value(t) - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
        forAll(cells, i)
        {
            const label cellI = cells[i];
            Udiag[cellI] +=
                V[cells[i]]*rho[cellI]*C0*C1m1b2*
                pow(magSqr(UpRel[cellI]), C1m1b2 - 1);
        }
    }
}

// ************************************************************************* //
