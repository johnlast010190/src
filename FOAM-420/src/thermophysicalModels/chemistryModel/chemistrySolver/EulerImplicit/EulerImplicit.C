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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistrySolver/EulerImplicit/EulerImplicit.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "matrices/simpleMatrix/simpleMatrix.H"

#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::EulerImplicit
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    chemistrySolver<ChemistryModel>(obr, phaseName),
    coeffsDict_(this->subDict("EulerImplicitCoeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter")),
    cTp_(this->nEqns()),
    obr_(obr),
    T_(const_cast<volScalarField&>(this->thermo().T())),
    p_(const_cast<volScalarField&>(this->thermo().p())),
    he_(const_cast<volScalarField&>(this->thermo().he())),
    mixCoeffs_(this->nSpecie_, 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::~EulerImplicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::EulerImplicit<ChemistryModel>::updateRRInReactionI
(
    const label index,
    const scalar pr,
    const scalar pf,
    const scalar corr,
    const label lRef,
    const label rRef,
    const scalar p,
    const scalar T,
    simpleMatrix<scalar>& RR
) const
{
    const Reaction<typename ChemistryModel::thermoType>& R =
        this->reactions_[index];

    forAll(R.lhs(), s)
    {
        const label si = R.lhs()[s].index;
        const scalar sl = R.lhs()[s].stoichCoeff;
        RR[si][rRef] -= sl*pr*corr;
        RR[si][lRef] += sl*pf*corr;
    }

    forAll(R.rhs(), s)
    {
        const label si = R.rhs()[s].index;
        const scalar sr = R.rhs()[s].stoichCoeff;
        RR[si][lRef] -= sr*pf*corr;
        RR[si][rRef] += sr*pr*corr;
    }
}


template<class ChemistryModel>
void Foam::EulerImplicit<ChemistryModel>::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT,
    label celli
) const
{
    const label nSpecie = this->nSpecie();
    simpleMatrix<scalar> RR(nSpecie, 0, 0);

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0, c[i]);
    }

    // Calculate the absolute enthalpy
    const scalar cTot = sum(c);

    //--------------------------------------

    scalar ha = 0.0;
    const word thermoType(this->specieThermo_[0].typeName());
    if (this->isMatModels_)
    {
        const scalar pOld = p_[celli];
        const scalar TOld = T_[celli];
        p_[celli] = p;
        T_[celli] = T;

        forAll(mixCoeffs_, specieI)
        {
            mixCoeffs_[specieI] = this->Y_[specieI][celli];
            this->Y_[specieI][celli] = (this->WModels_[specieI][0]*c[specieI]);
        }

        ha = this->HaMixtureModel_->operator[](celli);

        forAll(mixCoeffs_, specieI)
        {
            this->Y_[specieI][celli] = mixCoeffs_[specieI];
        }
        p_[celli] = pOld;
        T_[celli] = TOld;
    }
    else
    {
        typename ChemistryModel::thermoType mixture
        (
            (this->specieThermo_[0].W()*c[0])*this->specieThermo_[0]
        );
        for (label i=1; i<nSpecie; i++)
        {
            mixture += (this->specieThermo_[i].W()*c[i])*this->specieThermo_[i];
        }
        ha = mixture.Ha(p, T);
    }

    const scalar deltaTEst = min(deltaT, subDeltaT);

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pr, cr;
        label lRef, rRef;

        const scalar omegai =
            this->omegaI(i, c, T, p, pf, cf, lRef, pr, cr, rRef);

        scalar corr = 1;
        if (eqRateLimiter_)
        {
            if (omegai < 0)
            {
                corr = 1/(1 + pr*deltaTEst);
            }
            else
            {
                corr = 1/(1 + pf*deltaTEst);
            }
        }

        updateRRInReactionI(i, pr, pf, corr, lRef, rRef, p, T, RR);
    }

    // Calculate the stable/accurate time-step
    scalar tMin = GREAT;

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = 0;
        for (label j=0; j<nSpecie; j++)
        {
            d -= RR(i, j)*c[j];
        }

        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            const scalar cm = max(cTot - c[i], 1e-5);
            tMin = min(tMin, cm/d);
        }
    }

    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);

    // Add the diagonal and source contributions from the time-derivative
    for (label i=0; i<nSpecie; i++)
    {
        RR(i, i) += 1/deltaT;
        RR.source()[i] = c[i]/deltaT;
    }

    // Solve for the new composition
    c = RR.LUsolve();

    // Limit the composition
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0, c[i]);
    }


    //---------------------------------------------------------------------
    // Update the temperature
    if (this->isMatModels_)
    {
        forAll(mixCoeffs_, specieI)
        {
            mixCoeffs_[specieI] = this->Y_[specieI][celli];
            this->Y_[specieI][celli] = (this->WModels_[specieI][0]*c[specieI]);
        }
        const scalar pOld = p_[celli];
        const scalar TOld = T_[celli];
        const scalar HEOld = he_[celli];
        p_[celli] = p;
        T_[celli] = T;
        he_[celli] = ha;

        T = this->THaMixtureModel_->operator[](celli);

        forAll(mixCoeffs_, specieI)
        {
            this->Y_[specieI][celli] = mixCoeffs_[specieI];
        }
        p_[celli] = pOld;
        T_[celli] = TOld;
        he_[celli] = HEOld;
    }
    else
    {
        typename ChemistryModel::thermoType mixture
        (
            (this->specieThermo_[0].W()*c[0])*this->specieThermo_[0]
        );
        for (label i=1; i<nSpecie; i++)
        {
            mixture += (this->specieThermo_[i].W()*c[i])*this->specieThermo_[i];
        }
        T = mixture.THa(ha, p, T);
    }
}


// ************************************************************************* //
