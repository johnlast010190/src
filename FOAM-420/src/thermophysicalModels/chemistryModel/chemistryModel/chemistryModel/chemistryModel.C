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

#include "chemistryModel/chemistryModel/chemistryModel.H"
#include "mixtures/reactingMixture/reactingMixture.H"
#include "fields/Fields/uniformField/UniformField.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "include/dummyThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryModel<CompType, ThermoType>::chemistryModel
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    CompType(obr, phaseName),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_(CompType::template lookupOrDefault<scalar>("Treact", 0.0)),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_),
    T_(const_cast<volScalarField&>(this->thermo().T())),
    p_(const_cast<volScalarField&>(this->thermo().p())),
    isMatModels_(false)
{
    if (isType<dummyThermo>(specieThermo_[0]))
    {
        isMatModels_ = true;
        const objectRegistry& matObj = obr.subRegistry("materialModels");
        materialTables& mat = matObj.lookupObjectRef<materialTables>("materialTables");
        HcModels_.resize(nSpecie_);
        CpModels_.resize(nSpecie_);
        WModels_.resize(nSpecie_);
        HaModels_.resize(nSpecie_);
        forAll(HcModels_, i)
        {
            const word specieName(Y_[i].member());
            const word phaseName(Y_[i].group());
            HcModels_.set
            (
                i,
                mat.sTable(phaseName, specieName)[HcModel::typeName]
            );
            CpModels_.set
            (
                i,
                mat.sTable(phaseName, specieName)[CpModel::typeName]
            );
            WModels_.set
            (
                i,
                mat.sTable(phaseName, specieName)[WModel::typeName]
            );
            HaModels_.set
            (
                i,
                mat.sTable(phaseName, specieName)[HaModel::typeName]
            );
        }
        HaMixtureModel_ =
                mat.sTable(phaseName, word::null)[HaModel::typeName];

        THaMixtureModel_ =
                mat.sTable(phaseName, word::null)[THaModel::typeName];
    }

    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    obr.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }

    Info<< "chemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryModel<CompType, ThermoType>::~chemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    dcdt = Zero;

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            dcdt[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            dcdt[si] += sr*omegai;
        }
    }
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const Reaction<ThermoType>& R = reactions_[index];
    scalar w = omega(R, c, T, p, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c);
    const scalar kr = R.kr(kf, p, T, c);

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    for (label i = 0; i < nSpecie_; i++)
    {
        c_[i] = max(0.0, c[i]);
    }

    omega(c_, T, p, dcdt);

    // Constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
//    scalar cSum = 0.0;
    scalar dT = 0.0;
    scalar cp = 0.0;

    if (isMatModels_)
    {
        const scalar Told = T_[0];
        const scalar pold = p_[0];
        T_[0] = T;
        p_[0] = p;
        for (label i = 0; i < nSpecie_; i++)
        {
            const scalar W = WModels_[i][0];
//            cSum += c_[i];
            rho += W*c_[i];
        }
        for (label i=0; i<nSpecie_; i++)
        {
            cp += c_[i]*CpModels_[i][0]*WModels_[i][0];
        }
        cp /= rho;
        for (label i = 0; i < nSpecie_; i++)
        {
            scalar hi = HaModels_[i][0]*WModels_[i][0];
            dT += hi*dcdt[i];
        }
        T_[0] = Told;
        p_[0] = pold;
    }
    else
    {
        for (label i = 0; i < nSpecie_; i++)
        {
            const scalar W = specieThermo_[i].W();
//            cSum += c_[i];
            rho += W*c_[i];
        }
        for (label i=0; i<nSpecie_; i++)
        {
            cp += c_[i]*specieThermo_[i].cp(p, T);
        }
        cp /= rho;
        for (label i = 0; i < nSpecie_; i++)
        {
            scalar hi = specieThermo_[i].ha(p, T);
            dT += hi*dcdt[i];
        }
    }

    dT /= rho*cp;

    dcdt[nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0.0);
    }

    dfdc = Zero;

    // Length of the first argument must be nSpecie_
    omega(c_, T, p, dcdt);

    forAll(reactions_, ri)
    {
        const Reaction<ThermoType>& R = reactions_[ri];

        const scalar kf0 = R.kf(p, T, c_);
        const scalar kr0 = R.kr(kf0, p, T, c_);

        forAll(R.lhs(), j)
        {
            const label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c_[si] > SMALL)
                        {
                            kf *= el*pow(c_[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c_[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c_[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc(si, sj) -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc(si, sj) += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            const label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar er = R.rhs()[i].exponent;
                if (i == j)
                {
                    if (er < 1.0)
                    {
                        if (c_[si] > SMALL)
                        {
                            kr *= er*pow(c_[si] + VSMALL, er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c_[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= pow(c_[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc(si, sj) += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc(si, sj) -= sr*kr;
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1.0e-3;

    omega(c_, T + delta, p, dcdt_);
    for (label i=0; i<nSpecie_; i++)
    {
        dfdc(i, nSpecie_) = dcdt_[i];
    }

    omega(c_, T - delta, p, dcdt_);
    for (label i=0; i<nSpecie_; i++)
    {
        dfdc(i, nSpecie_) = 0.5*(dfdc(i, nSpecie_) - dcdt_[i])/delta;
    }

    dfdc(nSpecie_, nSpecie_) = 0;
    dfdc(nSpecie_ + 1, nSpecie_) = 0;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->obr(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    const label nReaction = reactions_.size();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T_[celli];
            const scalar pi = p_[celli];

            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] =
                    rhoi*Y_[i][celli]
                   /(isMatModels_ ? WModels_[i][0] : specieThermo_[i].W());
                cSum += c_[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];

                omega(R, c_, Ti, pi, pf, cf, lRef, pr, cr, rRef);

                forAll(R.rhs(), s)
                {
                    tc[celli] += R.rhs()[s].stoichCoeff*pf*cf;
                }
            }

            tc[celli] = nReaction*cSum/tc[celli];
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();
        forAll(Y_, i)
        {
            const scalar hi
            (
                isMatModels_ ? HcModels_[i][0] : specieThermo_[i].Hc()
            );
            forAll(Qdot, celli)
            {
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::chemistryModel<CompType, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    tmp<volScalarField::Internal> tRR
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    const Reaction<ThermoType>& R = reactions_[ri];

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T_[celli];
        const scalar pi = p_[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] =
                rhoi*Yi
               /(isMatModels_ ? WModels_[i][0] : specieThermo_[i].W());
        }

        const scalar omegai = omega
        (
            R,
            c_,
            Ti,
            pi,
            pf,
            cf,
            lRef,
            pr,
            cr,
            rRef
        );
        scalar dcdt = 0;
        forAll(R.lhs(), s)
        {
            if (R.lhs()[s].index == si)
            {
                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt -= sl*omegai;
            }
        }

        forAll(R.rhs(), s)
        {
            if (R.rhs()[s].index == si)
            {
                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt += sr*omegai;
            }
        }

        RR[celli] =
            dcdt*(isMatModels_ ? WModels_[si][0] : specieThermo_[si].W());
    }

    return tRR;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::chemistryModel<CompType, ThermoType>::calculateQdot
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tQdot = calculateRR(ri, si);
    tQdot->rename("Qdot");

    tQdot.ref() *=
        dimensionedScalar
        (
            "hi",
            dimEnergy/dimMass,
            isMatModels_ ? -HcModels_[si][0] : -specieThermo_[si].Hc()
        );

    return tQdot;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] =
                rhoi*Yi
               /(isMatModels_ ? WModels_[i][0] : specieThermo_[i].W());
        }

        omega(c_, Ti, pi, dcdt_);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] =
                dcdt_[i]
               *(isMatModels_ ? WModels_[i][0] : specieThermo_[i].W());
        }
    }
}


template<class CompType, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    CompType::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    scalarField c0(nSpecie_);

    // specify constant density value on demand
    scalar rho_0 = 0.0;
    if (this->thermo().found("rho_0"))
    {
        rho_0 = readScalar(this->thermo().lookup("rho_0", 0.0));
    }

    forAll(rho, celli)
    {
        scalar Ti = T_[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho_0 > 0 ? rho_0 : rho[celli];

            scalar pi = p_[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/(isMatModels_ ? WModels_[i][0] : specieThermo_[i].W());
                c0[i] = c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > SMALL)
            {
                scalar dt = timeLeft;
                this->solve(c_, Ti, pi, dt, this->deltaTChem_[celli], celli);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] =
                    (c_[i] - c0[i])
                   *(isMatModels_ ? WModels_[i][0] : specieThermo_[i].W())
                   /deltaT[celli];

                if (rho_0 > 0)
                {
                    RR_[i][celli] *= rho[celli]/rho_0;
                }
            }
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }
    }
    return deltaTMin;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
