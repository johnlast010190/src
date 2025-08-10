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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::fv::SolverOption

Description
    Finite volume solver template class to instantiate solverObjects as
    fvOptions.

SourceFiles
    SolverOption.C

\*---------------------------------------------------------------------------*/

#ifndef SolverOption_H
#define SolverOption_H

#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "solverObject/solverObject.H"
#include "solverOption.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class foamSolve;

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class SolverOption Declaration
\*---------------------------------------------------------------------------*/

template<class SolverObject>
class SolverOption
:
    public solverOption
{

protected:

    // Protected data

        SolverObject solverObject_;

        SolveTable<solveList> sourceDependencies_;
        HashTable<wordList> boundarySourceDependencies_;

public:

    //- Runtime type information
    TypeName("SolverOption");


    // Constructor

        //- Construct from components
        SolverOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

    // Destructor

        virtual ~SolverOption();

    // Member Functions

        // Access

            virtual class solverObject& solverObject() override
            {
                return solverObject_;
            }

        // Operations

            virtual void setSolveController(const foamSolve& solveController)
            override
            {
                solverObject_.setSolveController(solveController);
            }

            //- Initialise after all fvOptions created
            virtual bool initialise() override;

            //- Correct called from the relevant executeHook function
            // of fvOptionList
            virtual void correct() override;

            // Sources

                virtual void addSourceDependencies
                (
                    SolveTable<solveList>& dependencies
                ) override;

                virtual void addSup
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(eqn, fieldi);
                }

                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(eqn, fieldi);
                }

                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(eqn, fieldi);
                }

                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(eqn, fieldi);
                }

                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(eqn, fieldi);
                }

                virtual void addSup
                (
                    fvBlockMatrix<vector>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(eqn, fieldi);
                }

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(rho, eqn, fieldi);
                }

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(rho, eqn, fieldi);
                }

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(rho, eqn, fieldi);
                }

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(rho, eqn, fieldi);
                }

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(rho, eqn, fieldi);
                }

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvBlockMatrix<vector>& eqn,
                    const label fieldi
                ) override
                {
                    solverObject().addSup(rho, eqn, fieldi);
                }


                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const scalarField& pf,
                    scalarField& f,
                    scalarField& df
                ) override
                {
                    solverObject().addBoundarySource
                    (
                        fieldName, patchID, pf, f, df
                    );
                }

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const vectorField& pf,
                    vectorField& f,
                    vectorField& df
                ) override
                {
                    solverObject().addBoundarySource
                    (
                        fieldName, patchID, pf, f, df
                    );
                }

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const sphericalTensorField& pf,
                    sphericalTensorField& f,
                    sphericalTensorField& df
                ) override
                {
                    solverObject().addBoundarySource
                    (
                        fieldName, patchID, pf, f, df
                    );
                }

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const symmTensorField& pf,
                    symmTensorField& f,
                    symmTensorField& df
                ) override
                {
                    solverObject().addBoundarySource
                    (
                        fieldName, patchID, pf, f, df
                    );
                }

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const tensorField& pf,
                    tensorField& f,
                    tensorField& df
                ) override
                {
                    solverObject().addBoundarySource
                    (
                        fieldName, patchID, pf, f, df
                    );
                }


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFvSolverOption(SolverObject)                                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Foam::fv::SolverOption<Foam::fv::SolverObject>,                        \
        Foam::fv::SolverObject::typeName.c_str(),                              \
        0                                                                      \
    );                                                                         \
                                                                               \
    Foam::fv::option::adddictionaryConstructorToTable                          \
        <Foam::fv::SolverOption<Foam::fv::SolverObject>>                       \
        addSolverOption##SolverObject##dictionary##ConstructorTooptionTable_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolverOption.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
