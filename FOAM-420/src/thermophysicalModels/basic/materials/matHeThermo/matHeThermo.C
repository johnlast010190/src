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
    (c) 2021-2023 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matHeThermo.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"

#include "thermo/energy/matAbsoluteEnthalpy/matAbsoluteEnthalpy.H"
#include "thermo/energy/matAbsoluteInternalEnergy/matAbsoluteInternalEnergy.H"
#include "thermo/energy/matSensibleEnthalpy/matSensibleEnthalpy.H"
#include "thermo/energy/matSensibleInternalEnergy/matSensibleInternalEnergy.H"
#include "transport/alphahKappaOverCp/alphahKappaOverCp.H"
#include "transport/vAlphahKappaOverCp/vAlphahKappaOverCp.H"
#include "transport/strainRateLaw/strainRateLaw.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
void Foam::matHeThermo<BasicThermo, MixtureType>::heBoundaryCorrection
(
    fvPatchScalarField& pf
)
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
void Foam::matHeThermo<BasicThermo, MixtureType>::heBoundaryCorrection
(
    volScalarField& h
)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        heBoundaryCorrection(hBf[patchi]);
    }
}


template<class BasicThermo, class MixtureType>
Foam::materialTables&
Foam::matHeThermo<BasicThermo, MixtureType>::matLookupOrConstruct
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    if (!obr.foundObject<objectRegistry>("materialModels"))
    {
        obr.store(new materialTables(obr, dict));
    }
    else if
    (
        !obr.subRegistry("materialModels").foundObject<materialTables>
        (
            "materialTables"
        )
    )
    {
        obr.store(new materialTables(obr, dict));
    }

    return
        obr.subRegistry
        (
            "materialModels"
        ).lookupObjectRef<materialTables>("materialTables");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
void Foam::matHeThermo<BasicThermo, MixtureType>::reportFlags() const
{
    Info<< nl << "Solver settings from material properties:" << nl;

    // Isochoric?
    if (materials_.isochoric(this->phaseName_))
    {
        Info<< "\"isochoric\"";
    }
    else
    {
        Info<< "\"non-isochoric\"";
    }

    // Compressible?
    if (materials_.incompressible(this->phaseName_))
    {
        Info<< " \"incompressible\"";
    }
    else
    {
        Info<< " \"compressible\"";
    }

    // Constant Cpv?
    if (materials_.isCpvConst(this->phaseName_))
    {
        Info<< " \"constant ";
    }
    else
    {
        Info<< " \"varying ";
    }
    if (energyName(*this) == "h")
    {
        Info<< "Cp\"";
    }
    else if (energyName(*this) == "e")
    {
        Info<< "Cv\"";
    }
    else
    {
        Info<< "Cpv\"";
    }

    // Isotropic?
    if (materials_.isotropic(this->phaseName_))
    {
        Info<< " \"isotropic\"";
    }
    else
    {
        Info<< " \"non-isotropic\"";
    }

    // Buoyant?
    if (isBuoyant_)
    {
        if (isDistinctBuoyancy_)
        {
            Info<< " \"distinct buoyancy\"";
        }
        else
        {
            Info<< " \"buoyant\"";
        }
    }
    Info<< nl << endl;
}


template<class BasicThermo, class MixtureType>
void Foam::matHeThermo<BasicThermo, MixtureType>::initBuoyancyFlags
(
    const word& buoyantName
)
{
    matScalarTable& mat =
        materials_.sTable(this->phaseName_, word::null);
    if (mat.found("buoyancy"))
    {
        isBuoyant_ = true;
        if (mat["rho"] != mat["buoyancy"])
        {
            isDistinctBuoyancy_ = true;
        }
    }
}


template<class BasicThermo, class MixtureType>
void Foam::matHeThermo<BasicThermo, MixtureType>::init()
{
    materials_.addSpeciesAndSpeciesMixtures<specieModels>
    (
        this->phaseName_,
        materialsSpecie::typeName,
        thermodynamics::typeName
    );
    materials_.addSpeciesAndSpeciesMixtures<equationOfState>
    (
        this->phaseName_
    );
    materials_.addSpeciesAndSpeciesMixtures<thermodynamics>
    (
        this->phaseName_
    );
    materials_.addSpeciesAndSpeciesMixtures<energyConversion>
    (
        this->phaseName_,
        standardThermo::typeName,
        thermodynamics::typeName
    );

    // Add energy entry on the specie level (without dict entry)
    // (thermodynamics has to be available)
    materials_.addSpecies
    (
        energy::typeName,
        energy::listModels(),
        this->phaseName_,
        dictionary::lookup<word>("energy"),
        thermodynamics::typeName
    );

    // Add energy entry on the phase level (without dict entry)
    // (thermodynamics has to be available)
    if (this->phaseName_ != word::null)
    {
        materials_.addGroup
        (
            {energy::typeName},
            energy::listModels(),
            this->phaseName_,
            word::null,
            dictionary::lookup<word>("energy"),
            false,
            word::null,
            wordList::null()
        );
    }

    const dictionary& phaseDict =
        this->dictionary::optionalSubDict(this->phaseName_);
    if (word(phaseDict.lookup("materialType")) == "solid")
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            kappaModel::typeName + "Model",
            wordList({kappaModel::typeName, vKappaModel::typeName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({kappaModel::typeName, vKappaModel::typeName})
        );
        materials_.addSpeciesAndSpeciesMixtures
        (
            alphahModel::typeName + "Model",
            wordList({alphahModel::typeName}),
            this->phaseName_,
            alphahKappaOverCp::typeName,
            kappaModel::typeName + "Model",
            wordList({alphahModel::typeName})
        );
        materials_.addSpeciesAndSpeciesMixtures
        (
            vAlphahModel::typeName + "Model",
            wordList({vAlphahModel::typeName}),
            this->phaseName_,
            vAlphahKappaOverCp::typeName,
            kappaModel::typeName + "Model",
            wordList({vAlphahModel::typeName})
        );
    }
    else
    {
        if (he_.db().template foundObject<volVectorField>("U"))
        {
            // Strain rate model (requires U)
            materials_.addSpeciesAndSpeciesMixtures
            (
                strainRateModel::typeName + "Model",
                wordList({strainRateModel::typeName}),
                this->phaseName_,
                strainRateLaw::typeName,
                muModel::typeName + "Model",
                wordList({strainRateModel::typeName})
            );
        }
        else
        {
            Info<< "Material library didn't find velocity field."
                << " Non-newtonian models can't be used." << endl;
        }

        materials_.addSpeciesAndSpeciesMixtures
        (
            muModel::typeName + "Model",
            wordList({muModel::typeName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({muModel::typeName})
        );
        materials_.addSpeciesAndSpeciesMixtures
        (
            kappaModel::typeName + "Model",
            wordList({kappaModel::typeName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({kappaModel::typeName})
        );
        materials_.addSpeciesAndSpeciesMixtures
        (
            alphahModel::typeName + "Model",
            wordList({alphahModel::typeName}),
            this->phaseName_,
            alphahKappaOverCp::typeName,
            kappaModel::typeName + "Model",
            wordList({alphahModel::typeName})
        );
    }

    materials_.addSpeciesAndSpeciesMixtures<energy>(this->phaseName_);

    if (phaseDict.found(limitModel::typeName + "Model"))
    {
        materials_.addGroup
        (
            {limitModel::typeName + "Model"},
            {limitModel::typeName},
            this->phaseName_
        );
    }

    const word buoyantName("buoyancyModel");
    const bool buoyancyModelPresent
    (
        materials_.isModelInDict(buoyantName, this->phaseName_)
    );

    // Seach for distinct buoyant model only if the case is buoyant
    // and with buoyancy rho different from rhoModel
    if (buoyancyModelPresent)
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            buoyantName,
            wordList({"buoyancy"}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({rhoModel::typeName})
        );
    }

    materials_.addGroup
    (
        {TModel::typeName + "Model"},
        {
            TModel::typeName,
            THsModel::typeName,
            THaModel::typeName,
            TEsModel::typeName,
            TEaModel::typeName
        },
        this->phaseName_,
        word::null,
        energyInverseTemperature::typeName,
        false,
        word::null,
        wordList::null()
    );

    // Update model pointers and dependency lists
    materials_.linkModelsAndUpdateTables();

    // Linking models are added only on updateTables steps
    if (buoyancyModelPresent)
    {
        initBuoyancyFlags(buoyantName);
    }

    // Check for circular model dependencies
    materials_.checkDepenencies();

    if (this->phaseName_ == word::null)
    {
        // Output material information
        materials_.reportModelsInfo();

        // Report solver flags
        reportFlags();
    }

    // Initialize all the times of he
    initHe(this->p_, he_);
}


template<class BasicThermo, class MixtureType>
void Foam::matHeThermo<BasicThermo, MixtureType>::initHe
(
    const volScalarField& p,
    volScalarField& he
)
{
    he.forceAssign(materials_(HEModel::typeName, this->phaseName_)());
    this->heBoundaryCorrection(he);

    // Note: T does not have oldTime
    if (p.nOldTimes() > 0)
    {
        materials_.updateScalarField
        (
            this->phasePropertyName("p"),
            p.oldTime()
        );
        initHe(p.oldTime(), he.oldTime());
    }

    // Field always reset back to original
    if (this->p_.nOldTimes() > 0 && p.nOldTimes() == 0)
    {
        materials_.updateScalarField(this->phasePropertyName("p"), this->p_);
    }
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::tmpBoundaryField
(
    const scalarField& p,
    const scalarField& T,
    const label patchi,
    const word& modelName
) const
{
    scalarField Told, pOld;
    const bool isTsame = (&T == &this->T_.boundaryField()[patchi]);
    const bool isPsame = (&p == &this->p_.boundaryField()[patchi]);
    if (!isTsame)
    {
        Told = this->T_.boundaryField()[patchi];
        this->T_.boundaryFieldRef()[patchi].forceAssign(T);
    }
    if (!isPsame)
    {
        pOld = this->p_.boundaryField()[patchi];
        this->p_.boundaryFieldRef()[patchi].forceAssign(p);
    }

    tmp<scalarField> tField
    (
        materials_(modelName, this->phaseName_).boundaryField()[patchi]
    );

    // Reset T/p back to the original
    if (!isTsame)
    {
        this->T_.boundaryFieldRef()[patchi].forceAssign(Told);
    }
    if (!isPsame)
    {
        this->p_.boundaryFieldRef()[patchi].forceAssign(pOld);
    }

    return tField;
}


template<class BasicThermo, class MixtureType>
Foam::word Foam::matHeThermo<BasicThermo, MixtureType>::energyName
(
    const dictionary& dict
) const
{
    const word eName = dict.lookup<word>("energy");

    if (eName == matSensibleInternalEnergy::typeName)
    {
        return matSensibleInternalEnergy::name();
    }
    else if (eName == matSensibleEnthalpy::typeName)
    {
        return matSensibleEnthalpy::name();
    }
    else if (eName == matAbsoluteInternalEnergy::typeName)
    {
        return matAbsoluteInternalEnergy::name();
    }
    else if (eName == matAbsoluteEnthalpy::typeName)
    {
        return matAbsoluteEnthalpy::name();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown energy type: " << eName << nl
            << exit(FatalError);
    }

    return word::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::matHeThermo<BasicThermo, MixtureType>::matHeThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(obr, phaseName),
    MixtureType(*this, obr, phaseName),
    he_
    (
        IOobject
        (
            this->phasePropertyName(energyName(*this)),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        this->mesh(),
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    materials_(matLookupOrConstruct(obr, dict)),
    isBuoyant_(false),
    isDistinctBuoyancy_(false),
    solveTFromhe_(true)
{
    init();
}


template<class BasicThermo, class MixtureType>
Foam::matHeThermo<BasicThermo, MixtureType>::matHeThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermo(obr, phaseName),
    MixtureType(this->phaseDict(), obr, phaseName),
    he_
    (
        IOobject
        (
            this->phasePropertyName(energyName(*this)),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        this->mesh(),
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    materials_(matLookupOrConstruct(obr, *this)),
    isBuoyant_(false),
    isDistinctBuoyancy_(false),
    solveTFromhe_(true)
{
    init();
    BasicThermo::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::matHeThermo<BasicThermo, MixtureType>::~matHeThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
bool Foam::matHeThermo<BasicThermo, MixtureType>::incompressible() const
{
    return materials_.incompressible(this->phaseName_);
}


template<class BasicThermo, class MixtureType>
bool Foam::matHeThermo<BasicThermo, MixtureType>::isochoric() const
{
    return materials_.isochoric(this->phaseName_);
}


template<class BasicThermo, class MixtureType>
bool Foam::matHeThermo<BasicThermo, MixtureType>::isCpvConst() const
{
    return materials_.isCpvConst(this->phaseName_);
}


template<class BasicThermo, class MixtureType>
bool Foam::matHeThermo<BasicThermo, MixtureType>::isotropic() const
{
    return materials_.isotropic(this->phaseName_);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::matHeThermo<BasicThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    // Update field pointers in the models:
    materials_.updateScalarField(this->phasePropertyName("p"), p);
    materials_.updateScalarField(this->phasePropertyName("T"), T);

    tmp<volScalarField> tHe(materials_(HEModel::typeName, this->phaseName_)());

    // Reset model field pointers back
    materials_.updateScalarField(this->phasePropertyName("p"), this->p_);
    materials_.updateScalarField(this->phasePropertyName("T"), this->T_);

    return tHe;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    const scalarField pOld(this->p_, cells);
    const scalarField TOld(this->T_, cells);
    UIndirectList<scalar>(this->p_, cells) = p;
    UIndirectList<scalar>(this->T_, cells) = T;

    tmp<scalarField> tCellsHE(new scalarField(cells.size()));
    scalarField& cellsHE = tCellsHE.ref();
    const baseModels<scalar>& HeMatModel =
        materials_(HEModel::typeName, this->phaseName_);
    forAll(cells, i)
    {
        const label celli = cells[i];
        cellsHE[i] = HeMatModel[celli];
    }

    // Reset T/p back to the original
    UIndirectList<scalar>(this->T_, cells) = TOld;
    UIndirectList<scalar>(this->p_, cells) = pOld;

    return tCellsHE;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(p, T, patchi, HEModel::typeName);;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& he,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    const scalarField pOld(this->p_, cells);
    const scalarField TOld(this->T_, cells);
    const scalarField heOld(he_, cells);

    UIndirectList<scalar>(this->p_, cells) = p;
    UIndirectList<scalar>(this->T_, cells) = T0;
    UIndirectList<scalar>(he_, cells) = he;

    tmp<scalarField> tCellsT(new scalarField(cells.size()));
    scalarField& cellsT = tCellsT.ref();
    const baseModels<scalar>& TMatModel =
        materials_(TModel::typeName, this->phaseName_);
    forAll(cells, i)
    {
        const label celli = cells[i];
        cellsT[i] = TMatModel[celli];
    }

    // Reset T/p back to the original
    UIndirectList<scalar>(this->T_, cells) = TOld;
    UIndirectList<scalar>(this->p_, cells) = pOld;
    UIndirectList<scalar>(he_, cells) = heOld;

    return tCellsT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& he,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    scalarField Told, pOld, heOld;
    const bool isTsame = (&T0 == &this->T_.boundaryField()[patchi]);
    const bool isPsame = (&p == &this->p_.boundaryField()[patchi]);
    const bool isHEsame = (&he == &he_.boundaryField()[patchi]);
    if (!isTsame)
    {
        Told = this->T_.boundaryField()[patchi];
        this->T_.boundaryFieldRef()[patchi].forceAssign(T0);
    }
    if (!isPsame)
    {
        pOld = this->p_.boundaryField()[patchi];
        this->p_.boundaryFieldRef()[patchi].forceAssign(p);
    }
    if (!isHEsame)
    {
        heOld = he_.boundaryField()[patchi];
        const_cast<volScalarField&>(he_).boundaryFieldRef()[patchi].forceAssign(he);
    }

    tmp<scalarField> tField
    (
        materials_(TModel::typeName, this->phaseName_).boundaryField()[patchi]
    );

    // Reset T/p back to the original
    if (!isTsame)
    {
        this->T_.boundaryFieldRef()[patchi].forceAssign(Told);
    }
    if (!isPsame)
    {
        this->p_.boundaryFieldRef()[patchi].forceAssign(pOld);
    }
    if (!isHEsame)
    {
        const_cast<volScalarField&>(he_).boundaryFieldRef()[patchi].forceAssign(heOld);
    }

    return tField;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::hc() const
{
    return materials_(HcModel::typeName, this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(p, T, patchi, CpModel::typeName);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::Cp() const
{
    return materials_(CpModel::typeName, this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(p, T, patchi, CvModel::typeName);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::Cv() const
{
    return materials_(CvModel::typeName, this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(p, T, patchi, gammaModel::typeName);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::gamma() const
{
    return materials_(gammaModel::typeName, this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(p, T, patchi, CpvModel::typeName);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::Cpv() const
{
    return materials_(CpvModel::typeName, this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(p, T, patchi, CpByCpvModel::typeName);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::CpByCpv() const
{
    return materials_(CpByCpvModel::typeName, this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::kappa() const
{
    tmp<volScalarField> kappa(Cp()*this->alpha_);
    kappa.ref().rename("kappa");

    return kappa;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::matHeThermo<BasicThermo, MixtureType>::kappa
(
    const label patchi
) const
{
    return
        materials_(CpModel::typeName, this->phaseName_).boundaryField()[patchi]
       *this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<volScalarField> kappaEff(Cp()*(this->alpha_ + alphat));
    kappaEff.ref().rename("kappaEff");

    return kappaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        materials_(CpModel::typeName, this->phaseName_).boundaryField()[patchi]
       *(
           this->alpha_.boundaryField()[patchi]
         + alphat
        );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<volScalarField> talphaEff
    (
        materials_(CpByCpvModel::typeName, this->phaseName_)()
      *(this->alpha_ + alphat)
    );
    talphaEff.ref().rename("alphaEff");

    return talphaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        materials_
        (
            CpByCpvModel::typeName,
            this->phaseName_
        ).boundaryField()[patchi]
       *(
            this->alpha_.boundaryField()[patchi]
          + alphat
        );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::matHeThermo<BasicThermo, MixtureType>::buoyantRho() const
{
    return materials_("buoyancy", this->phaseName_)();
}


template<class BasicThermo, class MixtureType>
bool Foam::matHeThermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        materials_.updateSubDictsPtrs();
        return materials_.read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
