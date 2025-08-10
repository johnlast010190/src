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
    (c) 2017-2021 Esi Ltd.

Class
    Foam::fv::disperseEulerian

Description
    A solverObject for solving passive scalar transport.

SourceFiles
    disperseEulerian.C

\*---------------------------------------------------------------------------*/

#ifndef disperseEulerian_H
#define disperseEulerian_H

#include "solverObject/solverObject.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "cfdTools/general/include/fvCFD.H"
#include "sets/topoSets/cellSet.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "disperseEulerian/phase/phase.H"
#include "dynamicFvMesh/dynamicFvMesh.H"

#include "disperseEulerian/forceModels/dragModels/dragModel/dragModel.H"
#include "disperseEulerian/forceModels/liftModels/liftModel/liftModel.H"
#include "disperseEulerian/forceModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class disperseEulerian Declaration
\*---------------------------------------------------------------------------*/

class disperseEulerian
:
    public solverObject
{
protected:

     // Protected members
        bool solveTD_;
        bool alphaTD_;
        turbulenceModel* turb_;
        scalar stabDivAlpha_;

        // allow specifying U and phi
        word UName_;
        word phiName_;

        PtrList<decoupledEulerian::phase> phases_;
        PtrList<decoupledEulerian::dragModel> forcesD_;
        PtrList<decoupledEulerian::liftModel> forcesL_;
        PtrList<decoupledEulerian::turbulentDispersionModel> forcesTD_;

        // field stored [kg/m2] or [kg/m2 s]
        PtrList<volScalarField> sMass_;

        // field stored [kg] or [kg/s]
        PtrList<volScalarField> gMass_;

        Switch stochasticDispersionModel_;

        // register velocity flux for boundary conditions
        surfaceScalarField* phidPtr_;

        // Switch for time integration of flux in post processing
        Switch timeIntOn_;
        // initial time for integration in post processing
        scalar startTime_;
        // active post processing report
        bool reportOn_;

        label nNonOrthCorr_;

        //- Used to specify execution order in USF framework
        word fieldDependency_;

        //- Optional scaling factor for flow field
        scalar scalingFactor_;


public:

    //- Runtime type information
    TypeName("disperseEulerian");


    // Constructors

        //- Construct from components
        disperseEulerian
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    // Destructor

        virtual ~disperseEulerian();


    // Member Functions

            virtual void read(const dictionary& dict) override;

            virtual bool initialise() override;

            //- Identify any solves handled by this solverObject, their levels of
            // iteration and dependent solves
            virtual void getSolveGraph
            (
                wordList& solveNames,
                HashTable<wordList>& requiredDependencies,
                HashTable<wordList>& optionalDependencies,
                HashTable<wordList>& correctorMembers
            ) override
            {
                solveNames.append("alphad");
                if (fieldDependency_ != "none")
                {
                    requiredDependencies.insert("alphad", {fieldDependency_});
                }
                optionalDependencies.insert("alphad", {"fvMesh"});

                correctorMembers.insert
                (
                    solverObject::outerCorrectorName, solveNames
                );
            }

            //- Solve the equation system
            virtual void correct(const word&, const word&) override;

    // Access Functions

            //- Access continuous phase velocity
            tmp<volVectorField> U() const
            {
                if (UName_ == "")
                {
                    return turb_->U()*scalingFactor_;
                }
                else
                {
                    return obr_.lookupObject<volVectorField>(UName_)*scalingFactor_;
                }
            }

            //- Access k
            tmp<volScalarField> k() const
            {
                dimensionedScalar ds("ds", turb_->k()().dimensions(), VSMALL);
                return max(turb_->k()*pow(scalingFactor_,2), ds);
            }

            //- Access epsilon
            tmp<volScalarField> epsilon() const
            {
                dimensionedScalar ds("ds", turb_->epsilon()().dimensions(), SMALL);
                return max(turb_->epsilon()*pow(scalingFactor_,3), ds);
            }

            //- Access turbulent viscosity
            tmp<volScalarField> nut() const
            {
                dimensionedScalar ds("ds", turb_->nut()().dimensions(), SMALL);
                return max(turb_->nut()*scalingFactor_, ds);
            }

    // Access Functions through phases_

            //- Access disperse phase volume fraction
            const volScalarField& alphad(const label idx) const
            {
                return phases_[idx].alphad();
            }
            volScalarField& alphad(const label idx)
            {
                return phases_[idx].alphad();
            }

            //- Access disperse phase velocity
            const volVectorField& Ud(const label idx) const
            {
                return phases_[idx].Ud();
            }
            volVectorField& Ud(const label idx)
            {
                return phases_[idx].Ud();
            }

            const surfaceScalarField& phid(const label idx)
            {
                return *phidPtr_;
            }

    //- Additional access functions for source terms
            const PtrList<decoupledEulerian::phase>& phases() const
            {
                return phases_;
            }
            const PtrList<decoupledEulerian::dragModel>& dragModel() const
            {
                return forcesD_;
            }
            const PtrList<decoupledEulerian::liftModel>& liftModel() const
            {
                return forcesL_;
            }
            const PtrList<decoupledEulerian::turbulentDispersionModel>& tdModel() const
            {
                return forcesTD_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
