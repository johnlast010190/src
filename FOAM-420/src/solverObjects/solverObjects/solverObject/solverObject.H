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
    Foam::solverObject

Description
    Finite volume solver object base class to implement modular solvers,
    instantiated as an fvOption or functionObject

SourceFiles
    solverObject.C

\*---------------------------------------------------------------------------*/

#ifndef solverObject_H
#define solverObject_H

#include "fvMesh/fvMesh.H"
#include "meshes/MeshObject/MeshObject.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "VectorN/primitives/vector4/vector4.H"
#include "VectorN/primitives/VectorN.H"
#include "VectorN/primitives/TensorN.H"
#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "cfdTools/general/fvOptions/solveID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class foamSolve;

/*---------------------------------------------------------------------------*\
                      Class solverObject Declaration
\*---------------------------------------------------------------------------*/

class solverObject
:
    public MeshObject<fvMesh, UpdateableMeshObject, solverObject>
{
    //- Check if the solver object is in inactive region
    bool lookupInactiveRegions();

public:
    // Static data
    static const word timeLoopCorrectorName;
    static const word outerCorrectorName;

protected:

    // Protected data

        const word name_;

        const objectRegistry& obr_;

        fvMesh& mesh_;

        const word regionName_;

        //- Reference the he dictionary
        dictionary dict_;

        //- Reference to the optional '...Coeffs' subdictionary
        dictionary coeffs_;

        //- The phase/group name
        word phaseName_;

        //- Current corrector
        word correctorName_;

        //- Which corrector number we are on for each involved corrector
        HashTable<label> correctorNumber_;

        //- Whether on final iteration for each involved corrector
        HashTable<bool> finalIter_;

        //- Time start
        scalar timeStart_;

        //- Time end
        scalar timeEnd_;

        //- Solution interval - only solve every n outer iterations
        //  (0: never solve)
        autoPtr<Function1<scalar>> solveInterval_;

        //- The number of outer-correctors that have passed since the object was
        //  last active
        mutable label nCorrectorsSinceLastActive_;

        //- Force solving initial iteration regardless of solveInterval
        bool solveInitial_;

        //- Force solving final iteration regardless of solveInterval
        bool solveFinal_;

        //- Whether the last active solve returned finished = true
        bool lastActiveSolveWasSatisfied_;

        //- The fvOptions list - always valid after initialise() but may be
        // empty
        fv::options* fvOptionsPtr_;

        //- False if running from foamSolve,
        // true if from a legacy solver (default)
        bool legacySolver_;

        //- Keep track of outer corrector for legacy solve
        label legacyOuterCorrector_;

        //- Allow only one execution per time step (deprecated; only applies if legacy)
        Switch oneExecPerIter_;
        label lastTimeIndex_;

        //- Verbose output for case debugging
        Switch verbose_;

        //- Is object in inactive region?
        bool inactiveRegion_;

        //- The main foamSolve class
        const foamSolve* solveController_;

        //- Is an adjoint solver object
        bool adjoint_;

    // Protected functions



        fvMesh& getMesh(const objectRegistry&) const;

        //- Append phase name if present
        virtual word addPhaseName(const word& memberName)
        {
            return IOobject::groupName(memberName, phaseName_);
        }

        bool isStatic() const;

public:
        //- Return whether to solve this iteration
        virtual bool active() const;

    //- Runtime type information
    TypeName("solverObject");

    // Constructors

        //- Construct from components
        solverObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

    //- Destructor
    virtual ~solverObject();


    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const objectRegistry& obr() const
            {
                return obr_;
            }

            const dictionary& dict() const
            {
                return dict_;
            }

            const dictionary& coeffs() const
            {
                return coeffs_;
            }

            fv::options& fvOptions()
            {
                if (!fvOptionsPtr_)
                {
                    // Looks up or creates the fvOptions
                    fvOptionsPtr_ = &fv::options::New(mesh(), obr_);
                }
                return *fvOptionsPtr_;
            }

            //- Return start time
            const scalar& timeStart() const
            {
                return timeStart_;
            }

            //- Return end time
            const scalar& timeEnd() const
            {
                return timeEnd_;
            }

            void setTimeStart(const scalar& timeStart);

            void setTimeEnd(const scalar& timeEnd);

            void setSolveInitial(const bool setFlag)
            {
                solveInitial_ = setFlag;
            }

            //- Set whether running from legacy solver
            void setLegacy(const bool legacy)
            {
                legacySolver_ = legacy;
            }

            //- Return solve interval
            const Function1<scalar>& solveInterval() const
            {
                return solveInterval_;
            }

            //- Reset solve interval
            void resetSolveInterval(const Function1<scalar>& interval)
            {
                solveInterval_.reset(interval.clone().ptr());
            }

            //- Access to the main foamSolve class
            const foamSolve& solveController() const
            {
                FOAM_ASSERT(solveController_)
                {
                    FatalErrorInFunction
                        << "Solve controller not set." << exit(FatalError);
                }
                return *solveController_;
            }

            const bool& adjoint() const
            {
                return adjoint_;
            }

            void setSolveController(const foamSolve& solveController)
            {
                solveController_ = &solveController;
            }

        // Interface for derived classes. The ...Base functions call the over-
        // loaded methods and are not intended to be overloaded themselves.

            //- (Re)read controls dictionary. Called after construction and
            //  if file changes on disk
            virtual void read(const dictionary& dict)
            {}
            virtual void readBase(const dictionary& dict);

            //- Called to perform further initialisation after all solverObjects
            // and fvOptions in all regions have been created. Returning false
            // will cancel and destroy the solver object
            virtual bool initialise()
            {
                return true;
            }
            virtual bool initialiseBase();

            // Identify any solves handled by this solverObject. These can be
            // any unique identifier of each of the solution processes;
            // optionally the name of a field, but not necessarily. For each
            // solve, optionally also record the names of any additional fields
            // that are derived from the solve. Also return:
            //  * Any dependent solves (solve names or fields) which must be
            //    completed first.
            //  * The correctors which this solver object is part of.
            //  * For each corrector, return the solves which should be part of
            //    it. If a corrector name is given, all solves that are part
            //    of that corrector will be added.
            // The supplied tables are guaranteed empty when called.
            // Only override one of the three getSolveGraph() signatures below.
            // The first allows for dependent fields in different regions to be
            // specified, as well as derived fields; the second allows for
            // derived fields only; the third allows neither.
            virtual void getSolveGraph
            (
                wordList& solveNames,
                HashTable<wordList>& derivedFields,
                HashTable<solveList>& requiredDependencies,
                HashTable<solveList>& optionalDependencies,
                HashTable<wordList>& correctorMembers
            );
            virtual void getSolveGraph
            (
                wordList& solveNames,
                HashTable<wordList>& derivedFields,
                HashTable<wordList>& requiredDependencies,
                HashTable<wordList>& optionalDependencies,
                HashTable<wordList>& correctorMembers
            )
            {
                this->getSolveGraph
                (
                    solveNames,
                    requiredDependencies,
                    optionalDependencies,
                    correctorMembers
                );
            }
            virtual void getSolveGraph
            (
                wordList& solveNames,
                HashTable<wordList>& requiredDependencies,
                HashTable<wordList>& optionalDependencies,
                HashTable<wordList>& correctorMembers
            )
            {}
            virtual void getSolveGraphBase
            (
                wordList& solveNames,
                HashTable<wordList>& derivedFields,
                SolveTable<solveList>& requiredDependencies,
                SolveTable<solveList>& optionalDependencies,
                HashTable<wordList>& correctorMembers
            );

            //- Return maximum acceptable time step
            virtual scalar getMaxTimeStep()
            {
                return GREAT;
            }

            //- Called at the top of the corrector loop to determine if
            // enough correctors have been run at this level (including
            // the one about to be run, which is always run - possibly the
            // final corrector). Correctors will be continued until
            // all solvers have returned true here
            virtual bool isFinalCorrector
            (
                const label corrector,
                const word& correctorName
            )
            {
                return (correctorName != timeLoopCorrectorName);
            }
            virtual bool isFinalCorrectorBase
            (
                const label corrector,
                const word& correctorName
            );

            //- Called at the top of a corrector loop
            virtual void beginIteration
            (
                const label corrector,
                const word& correctorName,
                const bool finalIter
            )
            {}
            virtual void beginIterationBase
            (
                const label corrector,
                const word& correctorName,
                const bool finalIter
            );

            //- Called to return a matrix. Can return nothing, in which
            // case nothing will be solved.
            // finalSolve is initialised with the finalIteration data from the
            // mesh. Overwrite it to override the 'final solve' selection passed
            // to the solver. dictName is initialised with the word::null and
            // can be overwritten to specify the solver name to be selected in
            // fvSolution.
            virtual tmp<fvScalarMatrix> assembleScalarMatrix
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            )
            {
                return tmp<fvScalarMatrix>();
            }
            virtual tmp<fvScalarMatrix> assembleScalarMatrixBase
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            );
            virtual tmp<fvVectorMatrix> assembleVectorMatrix
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            )
            {
                return tmp<fvVectorMatrix>();
            }
            virtual tmp<fvVectorMatrix> assembleVectorMatrixBase
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            );
            virtual tmp<fvSymmTensorMatrix> assembleSymmTensorMatrix
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            )
            {
                return tmp<fvSymmTensorMatrix>();
            }
            virtual tmp<fvSymmTensorMatrix> assembleSymmTensorMatrixBase
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            );
            virtual tmp<fvTensorMatrix> assembleTensorMatrix
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            )
            {
                return tmp<fvTensorMatrix>();
            }
            virtual tmp<fvTensorMatrix> assembleTensorMatrixBase
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            );

            virtual tmp<fvBlockMatrix<vector4>> assembleVector4Matrix
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            )
            {
                return tmp<fvBlockMatrix<vector4>>();
            }

            virtual tmp<fvBlockMatrix<vector4>> assembleVector4MatrixBase
            (
                const word& solveName,
                bool& finalSolve,
                word& dictName
            );

            //- Called after the solve, if any. If the object wishes to do the
            // solve itself, it can do it here
            virtual void correct
            (
                const word& solveName, const word& regionName
            )
            {};


            virtual void correctBase
            (
                const word& solveName, const word& regionName
            );


            //- Called at the bottom of the corrector loop
            virtual void endIteration
            (
                const label corrector,
                const word& correctorName,
                const bool finalIter
            )
            {}
            virtual void endIterationBase
            (
                const label corrector,
                const word& correctorName,
                const bool finalIter
            );


            //- Called from correct() of SolverOption (via the hookOp) or
            // execute() of SolverFunction. Will assemble, solve, and correct
            // the solverObject for legacy solvers, to at least give a
            // partitioned solve via the same interface
            void solve();


        // Sources

            //- Derived objects should return the list of other equations to
            //  which this solver object contributes terms via fvOptions.
            //  This must be the name of the field that the equation is
            //  attached to, in order to be picked up by fvOptions.
            //  Dependencies should be filled with the dependent solves for each
            //  of these
            //  The first option is for the current region; the second allows
            //  source terms to be contributed to equations in different regions
            virtual void getSourceGraph
            (
                wordList& fieldNames,
                HashTable<wordList>& sourceDependencies
            )
            {}
            virtual void getSourceGraph
            (
                solveList& fieldNames,
                SolveTable<solveList>& sourceDependencies
            )
            {
                // This default implementation just adds current region name to
                // the single-region version
                wordList fNames;
                HashTable<wordList> sourceDeps;
                getSourceGraph(fNames, sourceDeps);
                fieldNames.resize(fNames.size());
                forAll(fNames, fi)
                {
                    fieldNames[fi] = solveID(fNames[fi], regionName_);
                }
                forAllConstIters(sourceDeps, depIter)
                {
                    solveID newKey = solveID(depIter.key(), regionName_);
                    sourceDependencies.insert
                    (
                        newKey,
                        solveList(depIter().size())
                    );
                    solveList& depList = sourceDependencies[newKey];
                    forAll(depIter(), depi)
                    {
                        depList[depi] = solveID(depIter()[depi], regionName_);
                    }
                }
            }

            virtual void getSourceGraphBase
            (
                solveList& fieldNames,
                SolveTable<solveList>& sourceDependencies
            )
            {
                FOAM_ASSERT
                (
                    !fieldNames.size() && !sourceDependencies.size()
                )
                {
                    FatalErrorInFunction
                        << "Unexpected parameters received."
                        << exit(FatalError);
                }

                getSourceGraph(fieldNames, sourceDependencies);
            }


            //- Add terms to the right hand side of the relevant equation
            //  returned in sourceFields

                virtual void addSup
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    fvBlockMatrix<vector>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                )
                {}

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvBlockMatrix<vector>& eqn,
                    const label fieldi
                )
                {}


        // Boundary sources

            //- Derived objects should return the list of patch IDs for any
            //  fields (key of the hash table) for which this solverObject
            //  contributes a boundary source. These are only added via a
            //  boundary condition that supports it.
            //  Dependencies should be filled with the dependent solves for each
            //  of these fields.
            virtual void getBoundarySourceGraph
            (
                HashTable<labelList>& fieldPatchIDs,
                HashTable<wordList>& boundarySourceDependencies
            )
            {}

            //- Add boundary source and derivative (optional)

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchIdx,
                    const scalarField& pf,
                    scalarField& f,
                    scalarField& df
                )
                {}

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const vectorField& pf,
                    vectorField& f,
                    vectorField& df
                )
                {}

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const sphericalTensorField& pf,
                    sphericalTensorField& f,
                    sphericalTensorField& df
                )
                {}

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const symmTensorField& pf,
                    symmTensorField& f,
                    symmTensorField& df
                )
                {}

                virtual void addBoundarySource
                (
                    const word& fieldName,
                    const label patchID,
                    const tensorField& pf,
                    tensorField& f,
                    tensorField& df
                )
                {}


        // Mesh update notifications

            //- Override to be notified when mesh moved
            virtual bool movePoints() override
            {
                return true;
            };

            //- Override to be notified when mesh underwent topology change
            virtual void updateMesh(const mapPolyMesh& mpm) override
            {};


        // Writing

            //- Write data - called by the functionObject according to its
            // write interval
            virtual void write()
            {}

            //- Write data for the IOobject - not needed
            virtual bool writeData(Ostream& os) const override
            {
                return true;
            }

            virtual void updateSolverObjectTimes()
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
