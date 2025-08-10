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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulenceSolver.H"
#include "solverOption/SolverOption.H"
#include "fluidThermo/fluidThermo.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "regionProperties/regionProperties.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbulenceSolver, 0);
}
}

makeFvSolverOption(turbulenceSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::turbulenceSolver::turbulenceSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    turbulencePtr_(nullptr),
    laminarTransportPtr_(nullptr),
    solnControlPtr_(nullptr),
    singleRegion_(true)
{
    //-this is a hack. Currently the turbulence uses the correct function to
    // solve the equations. Thus we cant control where to print the messame for
    // multiregion or singleregion.
    // Thus I reread the regionProperties here to check.
    // MultiInstance will not work with this. Needs to have access to scheduler
    // inside the solverObjects
    {
        regionProperties rp(this->time());
        const hashedWordList& groupNames = rp.groupNames();
        label totalRegions = 0;
        forAll(groupNames, i)
        {
            totalRegions += rp[groupNames[i]].size();
        }
        if (totalRegions>1)
        {
            singleRegion_ = false;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::turbulenceSolver::initialise()
{
    // First check for existing turbulence model
    turbulencePtr_ =
        obr_.lookupObjectRefPtr<turbulenceModel>
        (
            IOobject::groupName(turbulenceModel::propertiesName, phaseName_)
        );
    if (turbulencePtr_)
    {
        incompressible::turbulenceModel* icoTurbPtr =
            dynamic_cast<incompressible::turbulenceModel*>(turbulencePtr_);
        if (icoTurbPtr)
        {
            laminarTransportPtr_ =
                const_cast<transportModel*>(&icoTurbPtr->transport());
        }
    }
    else
    {
        // Create and store. Only create compressible version here for USF
        fluidThermo* fThermoPtr =
            dynamic_cast<fluidThermo*>
            (
                &multiphaseThermo::lookupOrCreate(obr_, phaseName_)
            );
        if (fThermoPtr)
        {
            Info<< "Creating turbulence model\n" << endl;
            turbulencePtr_ =
            (
                compressible::turbulenceModel::New
                (
                    obr_.lookupObject<volScalarField>
                    (
                        IOobject::groupName("rho", phaseName_)
                    ),
                    obr_.lookupObject<volVectorField>
                    (
                        IOobject::groupName("U", phaseName_)
                    ),
                    obr_.lookupObject<surfaceScalarField>
                    (
                        IOobject::groupName("phi", phaseName_)
                    ),
                    *fThermoPtr
                ).ptr()
            );
            turbulencePtr_->store();
        }
        else
        {
            // Solid thermo
            return false;
        }
    }
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);
    turbulencePtr_->validate();
    return true;
}


void Foam::fv::turbulenceSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames = {"turbulence"};
    optionalDependencies.insert("turbulence", {"materialProperties"});
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::turbulenceSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    // Legacy checks of the turbCorr() member in solution controls.
    // Generic 'solveFinal' keyword in solver object settings should be used
    // in preference to these settings
    pimpleControl* pimpleCtrl = dynamic_cast<pimpleControl*>(solnControlPtr_);
    foamCoupledControl* coupledCtrl =
        dynamic_cast<foamCoupledControl*>(solnControlPtr_);

    if
    (
        finalIter_[solverObject::outerCorrectorName]
    ||  (
            (!pimpleCtrl || !pimpleCtrl->turbOnFinalIterOnly())
         && (!coupledCtrl || !coupledCtrl->turbOnFinalIterOnly())
        )
    )
    {
        if (laminarTransportPtr_)
        {
            laminarTransportPtr_->correct();
        }
        if (!singleRegion_)
        {
            Info<< "Region: " << regionName << ":" << endl;
        }
        turbulencePtr_->correct();
    }
}

// ************************************************************************* //
