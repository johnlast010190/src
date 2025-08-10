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
    (c) 2019-2024 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fluidEnergySolver.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "interpolation/surfaceInterpolation/schemes/midPoint/midPoint.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "finiteVolume/ddtSchemes/steadyStateDdtScheme/steadyStateDdtScheme.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "equationOfState/BoussinesqLaw/BoussinesqLaw.H"
#include "equationOfState/rhoConstLaw/rhoConstLaw.H"
#include "materials/matHeThermo/matHeThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fluidEnergySolver, 0);
}
}

makeFvSolverOption(fluidEnergySolver);

const Foam::Enum<Foam::fv::fluidEnergySolver::formulationType>
Foam::fv::fluidEnergySolver::formulationTypeNames_
(
    {
        {formulationType::totalEnergyEnergy, "totalEnergyEnergy"},
        {formulationType::totalEnergyTemperature, "totalEnergyTemperature"},
        {formulationType::boussinesqEnergy, "boussinesqEnergy"},
        {
            formulationType::boussinesqTemperature,
            "boussinesqTemperature"
        }
    }
);


// * * * * * * * * * * * Private member functions  * * * * * * * * * * * * * //

void Foam::fv::fluidEnergySolver::createK()
{
    if (!K_.valid())
    {
        word UName = addPhaseName("U");
        const volVectorField& U =
            obr_.lookupObject<volVectorField>(UName);

        // Correct kinetic energy when the MRF is present
        if (hasMRF_)
        {
            volVectorField Urel("UrelK", U);
            fvOptions().makeRelative(Urel);
            K_.set(new volScalarField("K", 0.5*magSqr(Urel)));
        }
        else
        {
            K_.set(new volScalarField("K", 0.5*magSqr(U)));
        }

        if (U.nOldTimes())
        {
            const volVectorField* Uold = &U.oldTime();
            volScalarField* Kold = &K_().oldTime();
            if (hasMRF_)
            {
                volVectorField Urel("UrelK", *Uold);
                fvOptions().makeRelative(Urel);
                K_.set(new volScalarField("K", 0.5*magSqr(Urel)));
            }
            else
            {
                Kold->forceAssign(0.5*magSqr(*Uold));
            }

            while (Uold->nOldTimes())
            {
                Uold = &Uold->oldTime();
                Kold = &Kold->oldTime();
                if (hasMRF_)
                {
                    volVectorField Urel("UrelK", *Uold);
                    fvOptions().makeRelative(Urel);
                    K_.set(new volScalarField("K", 0.5*magSqr(Urel)));
                }
                else
                {
                    Kold->forceAssign(0.5*magSqr(*Uold));
                }
            }
        }
    }
}


const Foam::volScalarField& Foam::fv::fluidEnergySolver::updateK
(
    const volVectorField& U
)
{
    // Correct kinetic energy when the MRF is present
    if (hasMRF_)
    {
        volVectorField Urel("UrelK", U);
        fvOptions().makeRelative(Urel);
        K_() = 0.5*magSqr(Urel);
    }
    else
    {
        K_() = 0.5*magSqr(U);
    }

    return K_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fluidEnergySolver::fluidEnergySolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    solnControlPtr_(nullptr),
    thermoPtr_(nullptr),
    hasMRF_(false),
    lowFroudeStabilisation_
    (
        dict.lookupOrDefault<Switch>("lowFroudeStabilisation", true)
    ),
    lowFroudeStabilisationCoeff_
    (
        dict.lookupOrDefault<scalar>("lowFroudeStabilisationCoeff", 0.05)
    ),
    formulation_(totalEnergyEnergy)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::fluidEnergySolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    thermoPtr_ = &multiphaseThermo::lookupOrCreate(obr(), phaseName_);

    if (dict_.found("formulation"))
    {
        formulation_ = formulationTypeNames_.lookup("formulation", dict_);
    }
    else
    {
        if (isA<multiphaseThermo>(*thermoPtr_))
        {
            formulation_ = totalEnergyTemperature;
        }
        else if
        (
            thermoPtr_->distinctBuoyancy()
         && thermoPtr_->he().member() == "e"
         && thermoPtr_->materials().sTable
            (
                word::null, word::null
            )["buoyancy"]->type() == BoussinesqLaw::typeName
         && thermoPtr_->materials().sTable
            (
                word::null, word::null
            )[rhoModel::typeName]->type() == rhoConstLaw::typeName
        )
        {
            formulation_ = boussinesqEnergy;
        }
        else
        {
            formulation_ = totalEnergyEnergy;
        }
    }

    if (formulation_ == totalEnergyEnergy)
    {
        Info<< "Solving for internal energy in total energy formulation"
            << endl;
        thermoPtr_->setSolveTFromhe();
    }
    else if (formulation_ == totalEnergyTemperature)
    {
        Info<< "Solving for temperature in total energy formulation" << endl;
        thermoPtr_->setSolveheFromT();
    }
    else if (formulation_ == boussinesqEnergy)
    {
        Info<< "Solving for internal energy in Boussinesq formulation"
            << endl;
        thermoPtr_->setSolveTFromhe();
    }
    else if (formulation_ == boussinesqTemperature)
    {
        Info<< "Solving for temperature in Boussinesq formulation" << endl;
        thermoPtr_->setSolveheFromT();
    }
    else
    {
        NotImplemented;
    }

    mesh_.schemes().setFluxRequired(thermoPtr_->he().name());
    mesh_.schemes().setFluxRequired(thermoPtr_->T().name());

    // Check if we have MRF
    hasMRF_ = false;
    for (int i=0; i < fvOptions().optionList::size(); ++i)
    {
        if (fvOptions().optionList::operator[](i).isMRF())
        {
            hasMRF_ = true;
        }
    }

    Info<< "Creating field K\n" << endl;
    createK();

    return true;
}


void Foam::fv::fluidEnergySolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    word heTName;
    if (thermoPtr_->solvesTFromhe())
    {
        heTName = thermoPtr_->he().name();
    }
    else
    {
        heTName = thermoPtr_->T().name();
    }
    solveNames.append({heTName});

    // Dependencies
    optionalDependencies.insert
    (
        heTName,
        {addPhaseName("U"), thermoPtr_->p().name(), "fvMesh"}
    );

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);

    optionalDependencies.insert("materialProperties", solveNames);
}


bool Foam::fv::fluidEnergySolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else
    {
        return false;
    }
}


tmp<fvScalarMatrix>
Foam::fv::fluidEnergySolver::assembleScalarMatrix
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    basicThermo& thermo = *thermoPtr_;
    dimensionedScalar pRef(thermo.pRef());
    if (formulation_ == totalEnergyEnergy)
    {
        // Conservative formulation written in terms of energy
        // Typically used in single-phase applications

        volScalarField& he = thermo.he();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));
        const volVectorField& U =
            obr_.lookupObject<volVectorField>(addPhaseName("U"));

        // Update turbulence kinetic energy
        const volScalarField& K = updateK(U);

        const volScalarField& p = thermo.p();
        autoPtr<volScalarField> dpdt;
        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, he)());
        if (!steady && thermoPtr_->dpdt())
        {
            // Update pressure time derivative if needed
            dpdt.set(fvc::ddt(p).ptr());

            if (mesh_.moving())
            {
                dpdt() -= fvc::div(fvc::meshPhi(rho, U), p);
            }
        }

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<fvScalarMatrix> EEqn
        (
            fvm::div(phi, he) + fvc::div(phi, K)
         ==
            fvOptions()(rho, he)
        );

        if (!steady)
        {
            EEqn.ref() += fvm::ddt(rho, he) + fvc::ddt(rho, K);
        }

        if (he.member() == "e")
        {
            EEqn.ref() +=
                fvc::div
                (
                    fvc::absolute(phi, rho, U), (p + pRef)/rho, "div(phiv,p)"
                );
        }
        else if (!steady && thermoPtr_->dpdt())
        {
            EEqn.ref() += -dpdt();
        }

        // Add the term laplacian(kappaEff, T) to RHS
        if (thermo.isCpvConst())
        {
            // Faster calculation
            tmp<volScalarField> alphaEff(turbulence.alphaEff());
            word scheme = "laplacian(kappa";
            if (alphaEff().group() != word::null)
            {
                scheme += "." + alphaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -= fvm::laplacian(alphaEff(), he, scheme);
        }
        else
        {
            tmp<volScalarField> kappaEff(turbulence.kappaEff());
            word scheme = "laplacian(kappa";
            if (kappaEff().group() != word::null)
            {
                scheme += "." + kappaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -=
                changeVariable
                (
                    fvm::laplacian(kappaEff(), thermo.T(), scheme),
                    1/thermo.Cpv(),
                    he
                );
        }

        if (thermo.buoyant())
        {
            const uniformDimensionedVectorField& g =
                obr_.lookupObject<uniformDimensionedVectorField>("g");
            tmp<volScalarField> bRho
            (
                thermo.distinctBuoyancy()
              ? thermo.buoyantRho()
              : tmp<volScalarField>(rho)
            );
            EEqn.ref() -= bRho()*(U & g);

            if (lowFroudeStabilisation_ && lowFroudeStabilisationCoeff_ > 0)
            {
                // Add stabilisation to the diagonal for the above term
                dimensionedScalar gCmptSum
                (
                    "gCmptSum", g.dimensions(), cmptSum(cmptMag(g.value()))
                );
                tmp<volScalarField> offDiagCoeff
                (
                    lowFroudeStabilisationCoeff_*bRho*gCmptSum
                );
                offDiagCoeff->dimensions() *= U.dimensions()/he.dimensions();

                if (debug)
                {
                    volScalarField D("D", U.db(), mesh_, dimless);
                    D.primitiveFieldRef() = EEqn().D()()/mesh_.V();
                    volScalarField stabCoeff("stabCoeff", U.db(), mesh_, dimless);
                    stabCoeff.primitiveFieldRef() = offDiagCoeff();
                    if (mesh_.time().writeTime())
                    {
                        D.write();
                        stabCoeff.write();
                    }
                }

                // offDiagCoeff is the diagonal dominance that would be required
                // for these coupled terms. Ensure diagonal is at least as large
                offDiagCoeff->primitiveFieldRef() -=
                    min
                    (
                        max(EEqn().D()()/mesh_.V(), scalar(0)),
                        offDiagCoeff().primitiveField()
                    );
                EEqn.ref() += fvm::SuSp(offDiagCoeff(), he);
                EEqn.ref() -= fvc::Sp(offDiagCoeff, he);
            }
        }

        if (thermoPtr_->found("porosity"))
        {
            EEqn.ref() *= readScalar(thermoPtr_->lookup("porosity"));
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == totalEnergyTemperature)
    {
        // Conservative formulation written in terms of temperature. Typically
        // used for multiphase to avoid large discontinuities in energy.

        volScalarField& T = thermo.T();
        volScalarField& rho =
            obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));
        const surfaceScalarField& phiv =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phiv"));

        const volVectorField& U =
            obr_.lookupObject<volVectorField>(addPhaseName("U"));

        // Update turbulence kinetic energy
        const volScalarField& K = updateK(U);

        const volScalarField& p = thermo.p();

        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, T)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<volScalarField> Cv;
        if (isA<multiphaseThermo>(*thermoPtr_))
        {
            // Perform harmonic averaging of Cv consistent with averaging
            // temperature rather than energy
            Cv = 1/refCast<multiphaseThermo>(*thermoPtr_).rCv();
        }
        else
        {
            Cv = thermoPtr_->Cv();
        }
        tmp<fvScalarMatrix> EEqn
        (
            Cv()*fvm::div(phi, T)
          + (
                fvc::div
                (
                    fvc::absolute(phiv, U), p+pRef, "div(phiv,p)"
                )
              + fvc::ddt(rho, K) + fvc::div(phi, K)
            )
         == fvOptions()(rho*Cv(), T)
        );
        if (!steady)
        {
            EEqn.ref() += Cv*fvm::ddt(rho, T);
        }

        // Add the term laplacian(kappaEff, T) to RHS
        tmp<volScalarField> kappaEff(turbulence.kappaEff());
        word scheme = "laplacian(kappa";
        if (kappaEff().group() != word::null)
        {
            scheme += "." + kappaEff().group();
        }
        scheme += "," + thermo.T().name() + ")";
        EEqn.ref() -=
            fvm::laplacian(kappaEff(), thermo.T(), scheme);

        if (thermoPtr_->found("porosity"))
        {
            EEqn.ref() *= readScalar(thermoPtr_->lookup("porosity"));
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == boussinesqEnergy)
    {
        // Not conservative, and visous dissipation terms omitted which are
        // typically negligible in Boussinesq approximation

        volScalarField& he = thermo.he();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));
        const volVectorField& U =
            obr_.lookupObject<volVectorField>(addPhaseName("U"));
        const volScalarField& p = thermo.p();

        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, he)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<fvScalarMatrix> EEqn
        (
            fvm::div(phi, he)
        );

        if (!steady)
        {
            EEqn.ref() += fvm::ddt(rho, he);
        }

        if (he.name() == "h")
        {
            EEqn.ref() -=
                fvc::div
                (
                    fvc::absolute(phi, rho, U),
                    (p + pRef)/rho,
                    "div(phiv,p)"
                );
            if (!steady && thermoPtr_->dpdt())
            {
                volScalarField dpdt(fvc::ddt(p));
                if (mesh_.moving())
                {
                    dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
                }
                EEqn.ref() -= dpdt;
            }
        }

        EEqn.ref() -= fvOptions()(rho, he);

        // Add the term laplacian(kappaEff, T) to RHS
        if (thermo.isCpvConst())
        {
            // Faster calculation
            tmp<volScalarField> alphaEff(turbulence.alphaEff());
            word scheme = "laplacian(kappa";
            if (alphaEff().group() != word::null)
            {
                scheme += "." + alphaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -= fvm::laplacian(alphaEff(), he, scheme);
        }
        else
        {
            tmp<volScalarField> kappaEff(turbulence.kappaEff());
            word scheme = "laplacian(kappa";
            if (kappaEff().group() != word::null)
            {
                scheme += "." + kappaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -=
                changeVariable
                (
                    fvm::laplacian(kappaEff(), thermo.T(), scheme),
                    1/thermo.Cpv(),
                    he
                );
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == boussinesqTemperature)
    {
        // Not conservative, and visous dissipation terms omitted which are
        // typically negligible in Boussinesq approximation

        volScalarField& T = thermo.T();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));

        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, T)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<fvScalarMatrix> EEqn
        (
            fvm::div(phi, T)
        );

        if (!steady)
        {
            EEqn.ref() += fvm::ddt(rho, T);
        }

        tmp<volScalarField> Cv = thermo.Cv();
        EEqn.ref() *= Cv();
        EEqn.ref() -= fvOptions()(rho*Cv, T);

        // Add the term laplacian(kappaEff, T) to RHS
        tmp<volScalarField> kappaEff(turbulence.kappaEff());
        word scheme = "laplacian(kappa";
        if (kappaEff().group() != word::null)
        {
            scheme += "." + kappaEff().group();
        }
        scheme += "," + thermo.T().name() + ")";
        EEqn.ref() -= fvm::laplacian(kappaEff(), thermo.T(), scheme);

        if (thermoPtr_->found("porosity"))
        {
            EEqn.ref() *= readScalar(thermoPtr_->lookup("porosity"));
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else
    {
        NotImplemented;
    }
}


void Foam::fv::fluidEnergySolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    if
    (
        formulation_ == totalEnergyEnergy
     || formulation_ == boussinesqEnergy
    )
    {
        // Solved for energy

        volScalarField& he = thermoPtr_->he();

        if (verbose_)
        {
            Info<< he.name() << " min/max: "
                << min(he).value() << " "
                << max(he).value() << endl;
        }

        fvOptions().correct(he);
    }
    else if
    (
        formulation_ == totalEnergyTemperature
     || formulation_ == boussinesqTemperature
    )
    {
        // Solved for temperature

        volScalarField& he = thermoPtr_->he();
        volScalarField& T = thermoPtr_->T();

        fvOptions().correct(T);

        he.forceAssign(thermoPtr_->he(thermoPtr_->p(), T));

        // In multiphase case, set energies of phase thermos such as to produce
        // the same temperature for all
        if (isA<multiphaseThermo>(*thermoPtr_))
        {
            multiphaseThermo& mpThermo = refCast<multiphaseThermo>(*thermoPtr_);
            forAll(mpThermo.thermos(), ti)
            {
                basicThermo& thermo = mpThermo.thermos()[ti];
                thermo.he().forceAssign(thermo.he(thermo.p(), T));
            }
        }
    }
    else
    {
        NotImplemented;
    }
}


// ************************************************************************* //
