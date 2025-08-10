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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2010-2023 Esi Ltd

Class
    Foam::fv::option

Description
    Finite volume options abstract base class.  Provides a base set of
    controls, e.g.:
    \verbatim
        type            scalarExplicitSource    // source type
        active          on;                     // on/off switch
    \endverbatim

Note
    On evaluation, source/sink options are to be added to the equation R.H.S.

SourceFiles
    fvOption.C
    fvOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef fvOption_H
#define fvOption_H

#include "fvMatrices/fvMatricesFwd.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "VectorN/primitives/vector4/vector4.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "cfdTools/general/fvOptions/solveID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace fv
{
    class option;
    class optionList;

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class option
{
public:

    // Public enumerations

        //- Enumeration for function execution hook types
        enum execHookType
        {
            ehtConstrain,    // Run from the constrain() fvOption call
            ehtCorrect,      // Run from the field-correction functions
            ehtOperator,     // Run from the source operator() call
            ehtOuterCorrect, // Run from the outer-corrector loop
            ehtSolve,        // Run from the final outer-corrector
            ehtNone
        };


protected:

    // Protected data

        //- Source name
        const word name_;

        //- Model type
        const word modelType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        const objectRegistry& obr_;

        //- Top level source dictionary
        dictionary dict_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;

        //- Source active flag
        Switch active_;

        //- Reference frame flag - always apply frame transformations
        Switch MRF_;

        //- Field names to apply source to - populated by derived models
        wordList fieldNames_;

        //- If allocated, match the region name to apply the source to as well
        wordList regionNames_;

        //- When to solve
        const execHookType execHook_;

        //- Applied flag list - corresponds to each fieldNames_ entry
        List<bool> applied_;

        //- Field and patch names to apply boundary source to - populated by
        // derived models
        HashTable<labelList> boundarySourcePatchIDs_;

        //- Applied flag list for boundary sources - corresponds to each
        //  boundarySourcePatchIDs_ entry
        HashTable<List<bool>> boundaryApplied_;

        //- Get Mesh from give object registry
        const fvMesh& getMesh(const objectRegistry&) const;


public:

    //- Runtime type information
    TypeName("option");


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            option,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const objectRegistry& obr
            ),
            (name, modelType, dict, obr)
        );

        //- Word list of solver hook type names
        static const NamedEnum<execHookType, 6> execHookTypeNames_;


    // Constructors

        //- Construct from components
        option
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Return clone
        autoPtr<option> clone() const
        {
            NotImplemented;
        }

        //- Return pointer to new fvOption object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh
            const objectRegistry& obr_;

            const word& name_;

        public:

            iNew
            (
                const objectRegistry& obr,
                const word& name
            )
            :
                obr_(obr),
                name_(name)
            {}

            autoPtr<option> operator()(Istream& is) const
            {
                //const word name(is);
                const dictionary dict(is);

                return autoPtr<option>
                (
                    option::New(name_, dict, obr_)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<option> New
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~option();


    // Member Functions

        //- Perform initialisation after all options created
        //  Returning false will cause the option to be destroyed
        virtual bool initialise()
        {
            return true;
        };

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return dictionary
            inline const dictionary& coeffs() const;

            //- Return const access to the source active flag
            inline bool active() const;

            //- Return const access to the source active flag
            inline bool MRF() const;

            //- Set the applied flag to true for field index fieldi
            inline void setApplied(const label fieldi);

            //- Set the boundary applied flag to true for the field/patch index
            inline void setBoundaryApplied
            (
                const word& fieldName, const label fieldPatchi
            );

            //- Return const access to the execution hook
            const execHookType& execHook() const
            {
                return execHook_;
            }

            //- Return const access to the execution hook name
            inline const word execHookName() const;


        // Edit

            //- Return access to the source active flag
            inline Switch& active();

            //- Return access to the source active flag
            inline Switch& MRF();

        // Checks

            //- Is the source active?
            virtual bool isActive();

            //- Is the source an MRF?
            virtual bool isMRF() const;

            //- Return index of field name if found in fieldNames list
            virtual label applyToField
            (
                const word& fieldName, const word& regionName
            ) const;

            //- Return index of patch name if found in patchNames list
            virtual label applyToBoundaryFieldAndPatch
            (
                const word& fieldName,
                const label patchID
            ) const;

            //- Check that the source has been applied
            virtual void checkApplied() const;

            //- Check that the boundary source has been applied
            virtual void checkBoundaryApplied() const;


        // Evaluation

            //- Interface used by solverObjects to return source and
            //  boundary source dependencies
            virtual void addSourceDependencies
            (
                SolveTable<solveList>& dependencies
            )
            {}

            // Explicit and implicit sources

                virtual void addSup
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvBlockMatrix<vector>& eqn,
                    const label fieldi
                );


            // Explicit and implicit sources for compressible equations

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvBlockMatrix<vector>& eqn,
                    const label fieldi
                );


            // Explicit and implicit sources for phase equations

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Boundary flux sources

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


            // Constraints

                virtual void constrain
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvBlockMatrix<vector>&,
                    const label
                ){};

                virtual void constrain
                (
                    fvBlockMatrix<vector4>&,
                    const label
                ){};


            // Correction

                virtual void correct(volScalarField& field);
                virtual void correct(volVectorField& field);
                virtual void correct(volSphericalTensorField& field);
                virtual void correct(volSymmTensorField& field);
                virtual void correct(volTensorField& field);

                //- Non-field-specific correct
                virtual void correct();


            // Frame manipulations (only available to MRF=true options)

                //- Report call to frame function by non-frame object
                void frameError(word) const;

                //- Add the frame acceleration to the matrix
                virtual void addAcceleration
                (
                    const volVectorField& U,
                    volVectorField& ddtU
                ) const;

                //- Add the frame acceleration contribution to the momentum
                virtual void addAcceleration
                (
                    fvVectorMatrix& UEqn,
                    bool rhs = false
                ) const;

                //- Add the frame acceleration contribution to the momentum
                virtual void addAcceleration
                (
                    const volScalarField& rho,
                    fvVectorMatrix& UEqn,
                    bool rhs = false
                ) const;

                //- Add the frame acceleration contribution to the momentum
                virtual void addAcceleration
                (
                    fvBlockMatrix<vector>& UEqn
                ) const;

                //- Add the frame acceleration contribution to the momentum
                virtual void addAcceleration
                (
                    const volScalarField& rho,
                    fvBlockMatrix<vector>& UEqn
                ) const;

                //- Make the cell centre velocity relative
                virtual void makeRelative(volVectorField& U) const;

                //- Make the given absolute flux relative
                virtual void makeRelative(surfaceScalarField& phi) const;

                //- Make the given absolute boundary flux relative
                virtual void makeRelative
                (
                    FieldField<fvsPatchField, scalar>& phi
                ) const;

                //- Make the given absolute patch flux relative
                virtual void makeRelative
                (
                    Field<scalar>& phi,
                    const label patchi
                ) const;

                //- Make the given absolute mass-flux relative
                virtual void makeRelative
                (
                    const surfaceScalarField& rho,
                    surfaceScalarField& phi
                ) const;

                //- Make the given relative flux absolute
                virtual void makeAbsolute(surfaceScalarField& phi) const;

                //- Make the given relative mass-flux absolute
                virtual void makeAbsolute
                (
                    const surfaceScalarField& rho,
                    surfaceScalarField& phi
                ) const;

                //- Make the cell centre velocity absolute
                virtual void makeAbsolute(volVectorField& U) const;

                virtual void zero
                (
                    surfaceScalarField& phi
                ) const;

        // IO

            //- Write the source header information
            virtual void writeHeader(Ostream&) const;

            //- Write the source footer information
            virtual void writeFooter(Ostream&) const;

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/fvOptions/fvOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
