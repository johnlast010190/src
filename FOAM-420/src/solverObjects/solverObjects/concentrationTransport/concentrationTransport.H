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
    (c) 2010-2019 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::concentrationTransport

Description
    A function object for solving passive scalar transport. Where applicable
    optional entries include their default values. Some optional entries have
    no or NULL defaults

    Example of the function object specification:
    @verbatim
    tracer
    {
        type            concentrationTransport;
        region          region0;

        fieldName       Trace;
        phiName         phi; //optional
        rhoName         rho; //optional
        solver          smoothSolver
        {
            smoother    GaussSeidel;
            tolerance   1e-6;
            relTol      0;
        };              //optional

        solverOptions //this dictionary should always be present even if it is empty
        {
            solver          "smoothSolver{smoother GaussSeidel; tolerance 1e-6; relTol 0;}"; //optional - this is default if valid entry is not available

            solveTime       true; //optional
            timeScheme      Euler; //optional if valid default or explcicit entry is already present in fvSchemes
            relaxationFactor  1; //optional :  <1 necessary for steady state , default 0.7

            solveConvection true; //optional
            convectionScheme Gauss limitedLinear 1; //optional - this is the default if explicit entry is not in fvSchemes

            solveDiffusion  true; //optional
            diffusionScheme Gauss linear limited 0.333; //optional if valid default or explicit entry is already present in fvSchemes

            boundValue      0.0; //optional - bounding off by default
        }

        regions //any cell topo source - optional
        (
            boxToCell
            {
                box (0 0 -1) (0.1461 0.292 1);

                fixedValue 2;   //fixes cell value - optional - no default
                fixedSource 3;  //Su source value - optional - no default
                linearSource 4; //SuSp source value - optional - no default
            }
        );

        fieldDefinition // only necessary if the field is not defined on-disk
        {
            type scalar;
            dimensions [ 0 0 0 0 0 0 0 ];
            internalField uniform 0;

            boundaryTypeDefaults
            {

                walls {type zeroGradient;}

                outlets {type inletOutlet; inletValue uniform 0; value uniform 0;}

                inlets {type inletOutlet; inletValue uniform 0; value uniform 0;}
            }
            boundaryConditions
            {
                ABC
                {
                    type fixedValue; value uniform 10;
                }
                startTrace
                {
                    type fixedValue; value uniform 300;
                }
            }
        }
    }
    @endverbatim

SourceFiles
    concentrationTransport.C

\*---------------------------------------------------------------------------*/

#ifndef concentrationTransport_H
#define concentrationTransport_H

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
                          Class concentrationTransport Declaration
\*---------------------------------------------------------------------------*/

class concentrationTransport
:
    public solverObject
{
protected:

    // Protected data

        word fieldName_;

        word phiName_;

        word rhoName_;

        //- laminar diffusion
        autoPtr<dimensionedScalar> D_;

        //- turbulent Schmidt number
        autoPtr<dimensionedScalar> Sct_;

        //solver settings
        dictionary solverOptions_;

        //- cell sets for sources
        PtrList<cellSet> sourceCellSets_;

        //- fix the field value in specified cells
        labelList fixedValueSets_;
        scalarList fixedCellValues_;

        //- SuSp source coefficients
        labelList SuSpSets_;
        scalarList SuSpCoeffs_;

        //- Su sources
        labelList SuSets_;
        scalarList SuValues_;

        //- gravitational direction
        vector gHat_;

        //- terminal velocity
        dimensionedScalar uTerminal_;

        //- tracer density - for multiphase to prevent liquid
        // tracer entering gas
        autoPtr<word> phaseNamePtr_;

        //- Modify deltaT to speed-up passive transport (optional)
        label solveSpeedup_;

        //- Allow sub-cycling with specified number of time steps
        label nSubCycles_;

        //- Stabilize phasic transport (optional)
        dimensionedScalar residualPhaseField_;

        //- Used to specify execution order in USF framework
        word fieldDependency_;


    // Private functions
        //- Modify fvSchemes and fvSolution dictionaries to suite solver
        void setFvDicts();


        //- Disallow default bitwise copy construct
        //concentrationTransport(const concentrationTransport&);

        //- Disallow default bitwise assignment
        //void operator=(const concentrationTransport&);


public:

    //- Runtime type information
    TypeName("concentrationTransport");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        concentrationTransport
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    // Destructor

        virtual ~concentrationTransport();


    // Member Functions

        virtual bool initialise() override;

         //- Returns diffusivity
        tmp<volScalarField> diffusivity();

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
            solveNames.append(fieldName_);
            if (fieldDependency_ != "none")
            {
                requiredDependencies.insert(fieldName_, {fieldDependency_});
            }
            optionalDependencies.insert(fieldName_, {"fvMesh"});
            correctorMembers.insert
            (
                solverObject::outerCorrectorName, solveNames
            );
        }

        //- Update solver
        virtual void correct(const word&, const word&) override;

        //- Sample and write
        virtual void write() override;

        //- Read the parameters
        virtual void read(const dictionary&) override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
