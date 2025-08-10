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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "heThermo/heThermo.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
void Foam::heThermo<BasicThermo, MixtureType>::
heBoundaryCorrection(fvPatchScalarField& pf)
{
    if (isA<gradientEnergyFvPatchScalarField>(pf))
    {
        refCast<gradientEnergyFvPatchScalarField>(pf).gradient()
            = pf.fvPatchField::snGrad();
    }
    else if (isA<mixedEnergyFvPatchScalarField>(pf))
    {
        refCast<mixedEnergyFvPatchScalarField>(pf).refGrad()
            = pf.fvPatchField::snGrad();
    }
    else if (isA<blendedEnergyFvPatchScalarField>(pf))
    {
        blendedEnergyFvPatchScalarField& bepf =
            refCast<blendedEnergyFvPatchScalarField>(pf);
        heBoundaryCorrection(bepf.boundaryOne());
        heBoundaryCorrection(bepf.boundaryTwo());
    }
}


template<class BasicThermo, class MixtureType>
void Foam::heThermo<BasicThermo, MixtureType>::
heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        heBoundaryCorrection(hBf[patchi]);
    }
}


template<class BasicThermo, class MixtureType>
void Foam::heThermo<BasicThermo, MixtureType>::init
(
    const volScalarField& p,
    const volScalarField& T,
    volScalarField& he
)
{
    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& pCells = p.primitiveField();
    const scalarField& TCells = T.primitiveField();

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE
            (
                pCells[celli] + this->pRefValue(),
                TCells[celli]
            );
    }

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        heBf[patchi].forceAssign(
            this->he
            (
                p.boundaryField()[patchi],
                T.boundaryField()[patchi],
                patchi
            )
        );
    }

    this->heBoundaryCorrection(he);

    // Note: T does not have oldTime
    if (p.nOldTimes() > 0)
    {
        init(p.oldTime(), T.oldTime(), he.oldTime());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::heThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermo(obr, phaseName),
    MixtureType(*this, obr, phaseName),

    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init(this->p_, this->T_, he_);
}


template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::heThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(obr, dict, phaseName),
    MixtureType(*this, obr, phaseName),

    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init(this->p_, this->T_, he_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::~heThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> the
    (
        new volScalarField
        (
            IOobject
            (
                "he",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            he_.dimensions()
        )
    );

    volScalarField& he = the.ref();
    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TCells = T;

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE
            (
                pCells[celli] + this->pRefValue(), TCells[celli]
            );
    }

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        scalarField& hep = heBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(hep, facei)
        {
            hep[facei] =
                this->patchFaceMixture(patchi, facei).HE
                (
                    pp[facei] + this->pRefValue(), Tp[facei]
                );
        }
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    forAll(T, celli)
    {
        he[celli] = this->cellMixture(cells[celli]).HE
        (
            p[celli] + this->pRefValue(), T[celli]
        );
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    forAll(T, facei)
    {
        he[facei] =
            this->patchFaceMixture(patchi, facei).HE
            (
                p[facei] + this->pRefValue(), T[facei]
            );
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::hc() const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            he_.dimensions()
        )
    );

    volScalarField& hcf = thc.ref();
    scalarField& hcCells = hcf.primitiveFieldRef();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    volScalarField::Boundary& hcfBf = hcf.boundaryFieldRef();

    forAll(hcfBf, patchi)
    {
        scalarField& hcp = hcfBf[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp.ref();

    forAll(T, facei)
    {
        cp[facei] =
            this->patchFaceMixture(patchi, facei).Cp
            (
                p[facei] + this->pRefValue(), T[facei]
            );
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cp() const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp.ref();

    forAll(this->T_, celli)
    {
        cp[celli] =
            this->cellMixture(celli).Cp
            (
                this->p_[celli] + this->pRefValue(), this->T_[celli]
            );
    }

    volScalarField::Boundary& cpBf = cp.boundaryFieldRef();

    forAll(cpBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cpBf[patchi];

        forAll(pT, facei)
        {
            pCp[facei] =
                this->patchFaceMixture(patchi, facei).Cp
                (
                    pp[facei] + this->pRefValue(), pT[facei]
                );
        }
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv.ref();

    forAll(T, facei)
    {
        cv[facei] =
            this->patchFaceMixture(patchi, facei).Cv
            (
                p[facei] + this->pRefValue(), T[facei]
            );
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cv() const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cv = tCv.ref();

    forAll(this->T_, celli)
    {
        cv[celli] =
            this->cellMixture(celli).Cv
            (
                this->p_[celli] + this->pRefValue(), this->T_[celli]
            );
    }

    volScalarField::Boundary& cvBf = cv.boundaryFieldRef();

    forAll(cvBf, patchi)
    {
        cvBf[patchi] = Cv
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        );
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tgamma(new scalarField(T.size()));
    scalarField& gamma = tgamma.ref();

    forAll(T, facei)
    {
        gamma[facei] =
            this->patchFaceMixture(patchi, facei).gamma
            (
                p[facei] + this->pRefValue(), T[facei]
            );
    }

    return tgamma;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::gamma() const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> tgamma
    (
        new volScalarField
        (
            IOobject
            (
                "gamma",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimless
        )
    );

    volScalarField& gamma = tgamma.ref();

    forAll(this->T_, celli)
    {
        gamma[celli] =
            this->cellMixture(celli).gamma
            (
                this->p_[celli] + this->pRefValue(), this->T_[celli]
            );
    }

    volScalarField::Boundary& gammaBf = gamma.boundaryFieldRef();

    forAll(gammaBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = gammaBf[patchi];

        forAll(pT, facei)
        {
            pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
            (
                pp[facei] + this->pRefValue(),
                pT[facei]
            );
        }
    }

    return tgamma;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpv(new scalarField(T.size()));
    scalarField& Cpv = tCpv.ref();

    forAll(T, facei)
    {
        Cpv[facei] =
            this->patchFaceMixture(patchi, facei).Cpv
            (
                p[facei] + this->pRefValue(), T[facei]
            );
    }

    return tCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cpv() const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> tCpv
    (
        new volScalarField
        (
            IOobject
            (
                "Cpv",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& Cpv = tCpv.ref();

    forAll(this->T_, celli)
    {
        Cpv[celli] =
            this->cellMixture(celli).Cpv
            (
                this->p_[celli] + this->pRefValue(), this->T_[celli]
            );
    }

    volScalarField::Boundary& CpvBf = Cpv.boundaryFieldRef();

    forAll(CpvBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCpv = CpvBf[patchi];

        forAll(pT, facei)
        {
            pCpv[facei] =
                this->patchFaceMixture(patchi, facei).Cpv
                (
                    pp[facei] + this->pRefValue(), pT[facei]
                );
        }
    }

    return tCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpByCpv(new scalarField(T.size()));
    scalarField& CpByCpv = tCpByCpv.ref();

    forAll(T, facei)
    {
        CpByCpv[facei] =
            this->patchFaceMixture(patchi, facei).CpByCpv
            (
                p[facei] + this->pRefValue(), T[facei]
            );
    }

    return tCpByCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::CpByCpv() const
{
    const objectRegistry& obr = this->T_.db();
    tmp<volScalarField> tCpByCpv
    (
        new volScalarField
        (
            IOobject
            (
                "CpByCpv",
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimless
        )
    );

    volScalarField& CpByCpv = tCpByCpv.ref();

    forAll(this->T_, celli)
    {
        CpByCpv[celli] = this->cellMixture(celli).CpByCpv
        (
            this->p_[celli] + this->pRefValue(),
            this->T_[celli]
        );
    }

    volScalarField::Boundary& CpByCpvBf =
        CpByCpv.boundaryFieldRef();

    forAll(CpByCpvBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCpByCpv = CpByCpvBf[patchi];

        forAll(pT, facei)
        {
            pCpByCpv[facei] = this->patchFaceMixture(patchi, facei).CpByCpv
            (
                pp[facei] + this->pRefValue(),
                pT[facei]
            );
        }
    }

    return tCpByCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    tmp<scalarField> tT(new scalarField(h.size()));
    scalarField& T = tT.ref();

    forAll(h, celli)
    {
        T[celli] =
            this->cellMixture(cells[celli]).THE
            (
                h[celli], p[celli] + this->pRefValue(), T0[celli]
            );
    }

    return tT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{

    tmp<scalarField> tT(new scalarField(h.size()));
    scalarField& T = tT.ref();

    forAll(h, facei)
    {
        T[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).THE(h[facei], p[facei] + this->pRefValue(), T0[facei]);
    }

    return tT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappa() const
{
    tmp<Foam::volScalarField> kappa(Cp()*this->alpha_);
    kappa.ref().rename("kappa");
    return kappa;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::kappa
(
    const label patchi
) const
{
    return
        Cp
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> kappaEff(Cp()*(this->alpha_ + alphat));
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        Cp
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(
           this->alpha_.boundaryField()[patchi]
         + alphat
        );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpv()*(this->alpha_ + alphat));
    alphaEff.ref().rename("alphaEff");
    return alphaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->CpByCpv
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(
            this->alpha_.boundaryField()[patchi]
          + alphat
        );
}


template<class BasicThermo, class MixtureType>
bool Foam::heThermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        MixtureType::read(this->mesh(), *this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
