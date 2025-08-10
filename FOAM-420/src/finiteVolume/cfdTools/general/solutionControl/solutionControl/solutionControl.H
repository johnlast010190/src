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
    (c) 2016-2021 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd

Class
    Foam::solutionControl

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef solutionControl_H
#define solutionControl_H

#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solutionControl Declaration
\*---------------------------------------------------------------------------*/

class solutionControl
:
    public regIOobject
{
public:

    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


    // Static member functions

        template<class Type>
        static void maxTypeResidual
        (
            const objectRegistry& obr,
            const entry& solverPerfDictEntry,
            scalar& firstRes,
            scalar& lastRes
        );

        template<class Type>
        static void maxTypeBlockResidual
        (
            const objectRegistry& obr,
            const entry& solverPerfDictEntry,
            scalar& firstRes,
            scalar& lastRes
        );

        static scalar maxResidual
        (
            const objectRegistry& obr,
            const entry& solverPerfDictEntry,
            scalar& lastRes
        );

        //- Return index of field in residualControl_ if present
        virtual label applyToField
        (
            const word& fieldName,
            const bool useRegEx = true
        ) const;



protected:

    // Protected data

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Reference to the mesh database
        const objectRegistry& obr_;

        //- List of residual data per field
        List<fieldData> residualControl_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;


        // Solution controls

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indicate to solve using transonic algorithm
            bool transonic_;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            bool consistent_;

            //- Flag for Choi's modified momentum interpolation method
            bool modifiedMomentumInterp_;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            bool frozenFlow_;

            //- Flag to indicate that an adjoint optimization loop is active
            bool adjointOptimization_;

        // Evolution

            //- Current corrector loop index
            label corr_;

            //- Current non-orthogonal corrector loop index
            label corrNonOrtho_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void readControls(const bool absTolOnly);

        //- Return true if equation variables should be stored
        virtual bool storeVars() const;

        //- Store previous iteration field for vol<Type>Fields
        template<class Type>
        void storePrevIter() const;

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            const bool check = true,
            const bool force = false
        );

        const fvMesh& getMesh(const objectRegistry& obr) const
        {
            if (isA<fvSolutionRegistry>(obr))
            {
                return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
            }
            return dynamic_cast<const fvMesh&>(obr);
        }

private:

        //- Disallow default bitwise copy construct
        solutionControl(const solutionControl&);

        //- Disallow default bitwise assignment
        void operator=(const solutionControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solutionControl");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            solutionControl,
            dictionary,
            (fvMesh& mesh, const objectRegistry& obr),
            (mesh, obr)
        );


    // Constructors

        //- Construct from registry
        solutionControl(fvMesh& mesh, const word& algorithmName);

        //- Construct from registry
        solutionControl(fvMesh& mesh, const objectRegistry& obr, const word& algorithmName);

    // Selectors

        //- Run-time selector
        static autoPtr<solutionControl> New
        (
            fvMesh& mesh, const objectRegistry& obr
        );

        static autoPtr<solutionControl> NewSIMPLE
        (
            fvMesh& mesh, const objectRegistry& obr
        );

        //- Look up or construct if not found
        static solutionControl& lookupOrCreate
        (
            fvMesh& mesh, const objectRegistry& obr
        );

        //- Look up or construct a specific derived type
        template<class Type, class ...Args>
        static Type& lookupOrCreateType
        (
            fvMesh& mesh, const objectRegistry& obr, Args&&... args
        );

    //- Destructor
    virtual ~solutionControl();


    // Member Functions

        // Access

            //- Return the solution dictionary
            inline const dictionary& dict() const;

            //- Current corrector loop index
            inline label corr() const;

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;


        // Solution control

            //- Maximum number of outer correctors per timestep
            virtual label nOuterCorr() const
            {
                return 1;
            }

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            inline bool consistent() const;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            inline bool frozenFlow() const;

            //- Flag for Choi's modified momentum interpolation method
            inline bool modifiedMomentumInterp() const;

            inline bool adjointOptimization() const;

            void setAdjointOptimization();

            const List<fieldData>& residualControl() const
            {
                return residualControl_;
            }

            List<fieldData>& residualControl()
            {
                return residualControl_;
            }

            scalar maxResidual
            (
                const entry& solverPerfDictEntry,
                scalar& lastRes
            ) const;


        // Evolution

            //- Read controls from fvSolution dictionary
            virtual void readControls();

            //- Store previous iteration fields
            virtual void storePrevIterFields() const;

            //- Main control loop
            virtual bool loop() = 0;

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();

            //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() = 0;

        // Writing

            //- regIOobject write data function
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/solutionControl/solutionControl/solutionControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/solutionControl/solutionControl/solutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
