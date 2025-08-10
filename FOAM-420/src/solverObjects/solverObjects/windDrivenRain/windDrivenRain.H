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
    (c) 2017-2019 Esi Ltd.

Class
    Foam::fv::windDrivenRain

Description
    A fvOption for solving passive scalar transport.

SourceFiles
    windDrivenRain.C

\*---------------------------------------------------------------------------*/

#ifndef windDrivenRain_H
#define windDrivenRain_H

#include "solverObject/solverObject.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "cfdTools/general/include/fvCFD.H"
#include "sets/topoSets/cellSet.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class windDrivenRain Declaration
\*---------------------------------------------------------------------------*/

class windDrivenRain
:
    public solverObject
{
protected:

        //- Read fields from createFields.H
        dimensionedScalar Rh_;
        dimensionedScalar temp_;
        dimensionedScalar rhoa_;
        dimensionedScalar rhop_;
        dimensionedScalar mua_;
        scalar scalingFactor_;
        const scalar stabDivAlpha_;

        bool solveTD_;
        List<List<scalar>> phases_;

        volVectorField* UPtr_;
        //autoPtr<surfaceScalarField> phiwindPtr_;
        autoPtr<volScalarField> RePtr_;         // only needed to compute CdRe
        autoPtr<volScalarField> CdRePtr_;


        //- Read fields from createTDFields.H:
        PtrList<volScalarField> CtrainPtrL_;    // needed to compute nutrain and Rca
        autoPtr<volScalarField> kPtr_;          // needed to compute Ctrain and Rca
        autoPtr<volScalarField> epsilonPtr_;    // only needed to compute Ctrain
        autoPtr<volScalarField> nutPtr_;        // only needed to compute nutrain
        autoPtr<volScalarField> nutrainPtr_;    // needed in UEqns


        //- Read fields from createRainFields.H:
        PtrList<volVectorField> UrainPtrL_;
        PtrList<surfaceScalarField> phirainPtrL_;
        PtrList<volScalarField> alpharainPtrL_;
        PtrList<volScalarField> scrPtrL_;

        label nNonOrthCorr_;

        //- Used to specify execution order in USF framework
        word fieldDependency_;


    // Private functions

        //- compute material parameter
        void computeMaterialParam();

        //- compute drag coeff from Gunn & Kinzer table
        void computeCdRe();

public:

    //- Runtime type information
    TypeName("windDrivenRain");


    // Constructors

        //- Construct from components
        windDrivenRain
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    // Destructor

        virtual ~windDrivenRain();


    // Member Functions

            virtual void read(const dictionary& dict) override;

            //- Initialize solver (fields/variables)
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
                solveNames = {"alpharain"};
                if (fieldDependency_ != "none")
                {
                    requiredDependencies.insert("alpharain", {fieldDependency_});
                }
                optionalDependencies.insert("alpharain", {"fvMesh"});

                correctorMembers.insert
                (
                    solverObject::outerCorrectorName, solveNames
                );
            }

            //- Solve the equation system
            virtual void correct(const word&, const word&) override;


    // Access Functions

            //- Access velocity
            inline const volVectorField& U() const
            {
                if (!UPtr_)
                {
                    FatalErrorInFunction
                        << "U not allocated"
                        << nl << exit(FatalError);
                }
                return *UPtr_;
            }
            inline volVectorField& U()
            {
                if (!UPtr_)
                {
                    FatalErrorInFunction
                        << "U not allocated"
                        << nl << exit(FatalError);
                }
                return *UPtr_;
            }
/*
            //- Access wind flux
            inline const surfaceScalarField& phiwind() const
            {
                return phiwindPtr_();
            }
            inline surfaceScalarField& phiwind()
            {
                return phiwindPtr_();
            }
*/
            //- Access reynolds number
            inline const volScalarField& Re() const
            {
                return RePtr_();
            }
            inline volScalarField& Re()
            {
                return RePtr_();
            }

            //- Access drag times Reynolds
            inline const volScalarField& CdRe() const
            {
                return CdRePtr_();
            }
            inline volScalarField& CdRe()
            {
                return CdRePtr_();
            }

            //- Access rain Ct fields
            inline const volScalarField& Ctrain(const label& idx) const
            {
                return CtrainPtrL_[idx];
            }
            inline volScalarField& Ctrain(const label& idx)
            {
                return CtrainPtrL_[idx];
            }

            //- Access k
            inline const volScalarField& k() const
            {
                return kPtr_();
            }
            inline volScalarField& k()
            {
                return kPtr_();
            }

            //- Access epsilon
            inline const volScalarField& epsilon() const
            {
                return epsilonPtr_();
            }
            inline volScalarField& epsilon()
            {
                return epsilonPtr_();
            }

            //- Access turbulent viscosity
            inline const volScalarField& nut() const
            {
                return nutPtr_();
            }
            inline volScalarField& nut()
            {
                return nutPtr_();
            }

            //- Access turbulent viscosity of rain
            inline const volScalarField& nutrain() const
            {
                return nutrainPtr_();
            }
            inline volScalarField& nutrain()
            {
                return nutrainPtr_();
            }


            //- Access rain velocities
            inline const volVectorField& Urain(const label& idx) const
            {
                return UrainPtrL_[idx];
            }
            inline volVectorField& Urain(const label& idx)
            {
                return UrainPtrL_[idx];
            }

            //- Access rain fluxes
            inline const surfaceScalarField& phirain(const label& idx) const
            {
                return phirainPtrL_[idx];
            }
            inline surfaceScalarField& phirain(const label& idx)
            {
                return phirainPtrL_[idx];
            }

            //- Access rain phase fractions
            inline const volScalarField& alpharain(const label& idx) const
            {
                return alpharainPtrL_[idx];
            }
            inline volScalarField& alpharain(const label& idx)
            {
                return alpharainPtrL_[idx];
            }

            //- Access scr
            inline const volScalarField& scr(const label& idx) const
            {
                return scrPtrL_[idx];
            }
            inline volScalarField& scr(const label& idx)
            {
                return scrPtrL_[idx];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
