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
    (c) 2011-2016 OpenFOAM Foundation
    (c) Vuko Vukcevic, FMENA Zagreb.
    (c) Update by Hrvoje Jasak
    (c) 2019-2022 Esi Ltd.

Class
    Foam::fvBlockMatrix

Description
    fvBlockMatrix is an extension of fvMatrix for block coupled types.
    Support all Types (scalar, vector.. ) and VectorN
    Can be constructed/filled from:
    a) a segregated matrix fvMatrix
    b) another block matrix
    c) as a result of matrix operations (segregated/block)

SourceFiles
    fvBlockMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvBlockMatrix_H
#define fvBlockMatrix_H

#include "blockLduSystem/BlockLduSystem.H"
#include "fvMatrices/fvMatrices.H"
#include "fvMatrices/blockIndex/blockIndex.H"
#include "matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/blockLduSolvers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class fvBlockMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const fvBlockMatrix<Type>&);


/*---------------------------------------------------------------------------*\
                        Class fvBlockMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvBlockMatrix
:
    public BlockLduSystem<Type, Type>
{
    // Private data

        //- Internal field reference
        GeometricField<Type, fvPatchField, volMesh>& psi_;

        //- Class holding information about the indices
        blockIndex bIndex_;

    // Private member functions

        // Insertion functions for fvMatrix<Type2> into diagonal positions

            //- Insert internal field
            template<class fieldType>
            void insertSolutionVector
            (
                const direction dir,
                const Field<fieldType>& xSingle
            );

            //- Special treatment of BCs when we want components to behave
            // differently. E.g. pressureJump. (only 4th element of vector4
            // should have jump in AMI interpolation)
            template<class fieldType>
            void modifyBCCompSpecificMembers
            (
                const direction,
                const GeometricField<fieldType, fvPatchField, volMesh>&
            );

            //- Add matrix diagonal and source from fvMatrix
            template<class Type2>
            void insertDiagSource
            (
                const direction dir,
                fvMatrix<Type2>& matrix
            );

            //- Apply boundary manipulate or not
            //  Currently a hack. In Up-based systems only GIBs are using this
            //  structure.
            //  It will be used only when GIBs are present to avoif .ref or cast
            //  and having issues with cache
            template<class Type2>
            void applyBoundaryManipulate
            (
                fvMatrix<Type2>& matrix
            );

            //- Add matrix diagonal and source into this fvBlockMatrix
            template<class Type2>
            void insertDiagSource
            (
                const direction,
                const fvBlockMatrix<Type2>&
            );

            //- Add upper and lower part from fvMatrix
            template<class Type2>
            void insertUpperLower
            (
                const direction,
                const fvMatrix<Type2>&
            );

            //- Add upper and lower part from fvBlockMatrix
            template<class Type2>
            void insertUpperLower
            (
                const direction,
                fvBlockMatrix<Type2>&
            );

            //- Add coupling coeffs from fvMatrix
            template<class Type2>
            void updateCouplingCoeffs
            (
                const direction,
                const fvMatrix<Type2>&
            );

            //- Add coupling coeffs from fvBlockMatrix
            template<class Type2>
            void updateCouplingCoeffs
            (
                const direction,
                const fvBlockMatrix<Type2>&
            );

            //- Add boundary contributions from fvMatrix
            //  It is required if the algorithm later asks for flux (as SIMPLE)
            //  If so, the boundary contibutions must be stored in memory
            template<class Type2>
            void insertBoundaryCoeffsIfRequired
            (
                const direction,
                const fvMatrix<Type2>&
            );

            //- Add boundary contributions from fvBlockMatrix
            template<class Type2>
            void insertBoundaryCoeffsIfRequired
            (
                const direction,
                const fvBlockMatrix<Type2>&
            );

            //- Add boundary contributions from BlockLduSystem for the coupling
            //  matrices
            template<class blockType, class fieldType>
            void insertBoundaryCoeffsIfRequired
            (
                const direction,
                const direction,
                const BlockLduSystem<blockType, fieldType>&,
                const bool
            );

            //- Specialization for scalar fieldType
            //  Needed because generalised operator() for all the type
            // (for extracting components) does not exist in scalar type
            template<class blockType>
            void insertBoundaryCoeffsIfRequired
            (
                const direction,
                const direction,
                const BlockLduSystem<blockType, scalar>&,
                const bool
            );

            //- Store dimensions of the solving equations (not fields!)
            void checkAndSetDimensions
            (
                const direction,
                const direction,
                const dimensionSet&
            );

            //- Store dimensions of the solving equations (not fields!)
            void checkAndSetDimensions
            (
                const direction,
                const PtrList<dimensionSet>&
            );

            //- add the boundary contribution to the flux operation
            template<class Type2>
            void addBoundaryFlux
            (
                GeometricField<Type2, fvsPatchField, surfaceMesh>&,
                const direction
            ) const;

            //- set the snGrad corrections
            template<class Type2>
            void setFaceFluxCorrection
            (
                const GeometricField<Type2, fvsPatchField, surfaceMesh>&,
                const direction
            );

            //- construct the snGrad operator field
            void constructFaceFluxCorrection();

            //- add the snGrad corrections to the the flux operation
            template<class Type2>
            void addFaceFluxCorrection
            (
                GeometricField<Type2, fvsPatchField, surfaceMesh>&,
                const direction
            ) const;



        // Insertion functions for fvScalarMatrix into off-diagonal positions
        // (coupling matrices)

            //- Insert coupling matrix diag element into this fvBlockMatrix
            void insertCouplingDiag
            (
                const direction dirI,
                const direction dirJ,
                const scalarField& coeffIJ
            );

            //- Insert coupling matrix diag and source into this fvBlockMatrix
            void insertCouplingDiagSource
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );

            //- Insert coupling matrix lower and upper into this fvBlockMatrix
            void insertCouplingUpperLower
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );


        // Pressure - velocity (p-U) coupling specific functions

            //- Insert BlockLduSystem (obtained by implicit grad/div operator)
            // into this fvBlockMatrix
            template<class blockType, class fieldType>
            void insertBlock
            (
                const direction dirI,
                const direction dirJ,
                const BlockLduSystem<blockType, fieldType>& blockSystem,
                const bool incrementColumnDir
            );

            //- Insert source and coupling coeffs of BlockLduSystem
            // (eg. obtained by implicit grad/div operator)
            template<class blockType, class fieldType>
            void insertBoundaryContributions
            (
                const direction dirI,
                const direction dirJ,
                const BlockLduSystem<blockType, fieldType>& blockSystem,
                const bool incrementColumnDir
            );


public:

    // Constructors

        //- Construct given a field to solve for
        fvBlockMatrix(GeometricField<Type, fvPatchField, volMesh>&);

        //- Construct given a field to solve for and dimensions
        fvBlockMatrix
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const dimensionSet&
        );

        //- Construct as copy
        fvBlockMatrix(const fvBlockMatrix<Type>&);

        //- Construct as copy
        fvBlockMatrix(fvMatrix<Type>&);

        //- Construct as copy
        fvBlockMatrix(tmp<fvMatrix<Type>>);


    //- Destructor
    virtual ~fvBlockMatrix()
    {}


    // Member functions

        // Access

            //- Access to GeometricField
            GeometricField<Type, fvPatchField, volMesh>& psi()
            {
                return psi_;
            }

            //- Const reference to GeometricField
            const GeometricField<Type, fvPatchField, volMesh>& psi() const
            {
                return psi_;
            }

            //- Const reference to block index class
            const blockIndex& bIndex() const
            {
                return bIndex_;
            }

            //- Reference to block index class
            blockIndex& bIndex()
            {
                return bIndex_;
            }

            //- Function for fixing the solution in regions
            template<class Type2>
            void setValuesFromList
            (
                const word,
                const labelUList&,
                const Field<Type2>&
            );

            //- Specialization for scalars
            void setValuesFromList
            (
                const direction,
                const labelUList&,
                const scalarField&
            );

            //- Set reference for closed domains
            //  (like pressure bounded everywhere with Neumann-like BC)
            void setReference
            (
                const word,
                const label,
                const scalar&
            );


            //- Flux extraction for the Type (all equations)
            //  Careful! the BCs are not taken into consideration
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
                flux() const;


            //- Flux extraction for given field
            template<class Type2>
            tmp<GeometricField<Type2, fvsPatchField, surfaceMesh>> flux
            (
                const GeometricField<Type2, fvPatchField, volMesh>&
            ) const;


        // Insertion and retrieval public tools

            //- Retrieve part of internal field from this fvBlockMatrix
            template<class fieldType>
            void retrieveSolution
            (
                const direction,
                Field<fieldType>&
            ) const;

            //- Insert matrix into this fvBlockMatrix
            template<class Type2>
            void insertEquation
            (
                const direction&,
                fvMatrix<Type2>&
            );

            template<class Type2>
            void insertEquation
            (
                const direction&,
                tmp<fvMatrix<Type2>>
            );

            //- Insert fvBlockMatrix into this fvBlockMatrix
            template<class Type2>
            void insertEquation
            (
                const direction&,
                fvBlockMatrix<Type2>&
            );

            //- Insert existing block system (obtained by implicit grad/div
            // operator) into this fvBlockMatrix
            template<class blockType, class fieldType>
            void insertBlockCoupling
            (
                const direction dirI,
                const direction dirJ,
                const BlockLduSystem<blockType, fieldType>& blockSystem,
                const bool incrementColumnDir
            );

            //- Insert diagonal only equation coupling into this fvBlockMatrix
            //  Source compensation is done in function updateSourceCoupling()
            //  after all coupling terms are added.  HJ, 27/Apr/2015
            void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                const scalarField& coeffIJ
            );

            //- Insert scalar equation coupling into this fvBlockMatrix
            void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );

            //- Add field into block field
            void blockAdd
            (
                const direction dir,
                const scalarField& xSingle,
                Field<Type>& blockX
            );

            //- Update coupling of block system.
            //  Subtracts the block-coefficient coupling as specified by the
            //  user from the source, leaving the implicit update given by
            //  linearisation
            void updateSourceCoupling();

            //- Operators specific to NS-eq.
                //- Momentum - gradp BC couplings (like total-pressure BC)
                void addMomentumGradPCoupledBC(const volScalarField&);

        // Solver calls for fvBlockMatrix

            //- Solve returning the solution statistics.
            // Use the given solver controls
            BlockSolverPerformance<Type> solve(const dictionary&);

            //- Solve returning the solution statistics.
            // Solver controls read from fvSolution
            BlockSolverPerformance<Type> solve();

            //- Relax block equation
            void relax();

            //- Relax block equation
            void relax(const Field<Type>&);

            void relaxAndCheckDominance();

            //- Relax specific cells
            void relax(const scalar, const labelList&);

            //- Relax block equation
            void boundaryRelax();

            //- Relax block equation
            void boundaryProcRelax(const scalar);

            //- Relax block equation
            void boundaryRelax(const scalar, const label);

            //- Return the solver dictionary
            const dictionary& solverDict() const;

    // Member operators

        void operator=(const fvBlockMatrix<Type>&);
        void operator=(const tmp<fvBlockMatrix<Type>>&);

        void negate();

        void operator+=(const fvBlockMatrix<Type>&);
        void operator+=(const tmp<fvBlockMatrix<Type>>&);
        template<class Type2>
        void operator+=(fvBlockMatrix<Type2>&);
        template<class Type2>
        void operator+=(tmp<fvBlockMatrix<Type2>>&);

        void operator-=(const fvBlockMatrix<Type>&);
        void operator-=(tmp<fvBlockMatrix<Type>>&);
        template<class Type2>
        void operator-=(const fvBlockMatrix<Type2>&);
        template<class Type2>
        void operator-=(tmp<fvBlockMatrix<Type2>>&);

        void operator+=(const fvMatrix<Type>&);
        void operator+=(const tmp<fvMatrix<Type>>&);
        template<class Type2>
        void operator+=(const fvMatrix<Type2>&);
        template<class Type2>
        void operator+=(const tmp<fvMatrix<Type2>>&);

        void operator+=
        (
            const DimensionedField<Type, volMesh>&
        );

        void operator-=(const fvMatrix<Type>&);
        void operator-=(const tmp<fvMatrix<Type>>&);
        template<class Type2>
        void operator-=(const fvMatrix<Type2>&);
        template<class Type2>
        void operator-=(const tmp<fvMatrix<Type2>>&);

        void operator*=(const scalarField&);
        void operator*=(const tmp<scalarField>&);

        void operator*=(const scalar);

        // Operations used to fill the matrix using smaller sub-matrices

        template<class Type2>
        void operator()(const word, fvBlockMatrix<Type2>&);

        // Operations used to fill the coupling sub-matrices

        template<class Type2, class Type3>
        void operator()(const word, const word, BlockLduSystem<Type2, Type3>&);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvBlockMatrix<Type>&
        );
};

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const fvBlockMatrix<Type>&,
    const fvBlockMatrix<Type>&,
    const char*
);

template<class Type>
void checkMethod
(
    const fvBlockMatrix<Type>&,
    const fvMatrix<Type>&,
    const char*
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const fvBlockMatrix<Type>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const fvBlockMatrix<Type>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const dimensioned<Type>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const fvBlockMatrix<Type>&,
    const zero&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator==
(
    const tmp<fvBlockMatrix<Type>>&,
    const zero&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const DimensionedField<Type, volMesh>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const DimensionedField<Type, volMesh>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const fvBlockMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const tmp<fvBlockMatrix<Type>>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const dimensioned<Type>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator+
(
    const dimensioned<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const DimensionedField<Type, volMesh>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const DimensionedField<Type, volMesh>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const tmp<fvBlockMatrix<Type>>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const tmp<fvBlockMatrix<Type>>&
);


template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const fvBlockMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const tmp<fvBlockMatrix<Type>>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const dimensioned<Type>&,
    const fvBlockMatrix<Type>&
);

template<class Type>
tmp<fvBlockMatrix<Type>> operator-
(
    const dimensioned<Type>&,
    const tmp<fvBlockMatrix<Type>>&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.C"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrixOperations.C"
#endif

// Specialisation for scalars
#include "fvMatrices/fvBlockMatrix/fvScalarBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
