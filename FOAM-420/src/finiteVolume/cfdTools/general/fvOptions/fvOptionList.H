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
    (c) 2010-2023 Esi Ltd

Class
    Foam::fv::optionList

Description
    List of finite volume options

SourceFile
    optionList.C

\*---------------------------------------------------------------------------*/

#ifndef fvOptionList_H
#define fvOptionList_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

namespace fv
{
    class optionList;
}

Ostream& operator<<(Ostream& os, const fv::optionList& options);

namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class optionList Declaration
\*---------------------------------------------------------------------------*/

class optionList
:
    public PtrList<option>
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the region database
        const objectRegistry& obr_;

        //- Time index to check that all defined sources have been applied
        label checkTimeIndex_;

        //- optionLists for all regions
        UPtrList<optionList> regionOptions_;

        //- Names of all regions
        wordList regionNames_;


    // Protected Member Functions

        //- Return the "options" sub-dictionary if present otherwise return dict
        const dictionary& optionsDict(const dictionary& dict) const;

        //- Read options dictionary
        bool readOptions(const dictionary& dict);

        //- Check that all sources have been applied
        void checkApplied() const;

        //- Check that all boundary sources have been applied
        void checkBoundaryApplied() const;

        //- Get Mesh from give object registry
        const fvMesh& getMesh(const objectRegistry&) const;

        const word& thisRegionName() const
        {
            return obr_.name() == word::null ? mesh_.name() : obr_.name();
        }

        //- Disallow default bitwise copy construct
        optionList(const optionList&);

        //- Disallow default bitwise assignment
        void operator=(const optionList&);


public:

    //- Runtime type information
    TypeName("optionList");


    // Constructors

        //- Construct null
        optionList(const objectRegistry& obr);

        //- Construct from mesh and dictionary
        optionList(const objectRegistry& obr, const dictionary& dict);


    //- Destructor
    virtual ~optionList()
    {}


    // Member Functions

        //- Reset the source list
        void reset(const dictionary& dict, const objectRegistry& obj);

        //- Initialise after all fvOptions created
        virtual void initialise();

        //- Initialise and add an option (if initialisation returned success)
        bool addOption(autoPtr<option>& opt);

        // Sources

            //- Add the names of any dependent solves for sources
            void addSourceDependencies
            (
                SolveTable<solveList>& depList
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            template<class Type>
            tmp<fvBlockMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const bool block
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );

            template<class Type>
            tmp<fvBlockMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName,
                const bool block
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );


            template<class Type>
            tmp<fvBlockMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const bool block
            );

            template<class Type>
            tmp<fvBlockMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName,
                const bool block
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const geometricOneField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const geometricOneField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const geometricOneField& alpha,
                const geometricOneField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );


        // Boundary flux sources - return source and derivative

            template<class Type>
            tmp<Field<Type>> boundarySources
            (
                const word& fieldName,
                const label patchID,
                const Field<Type>& pf,
                Field<Type>& bsourceDerivs
            );


        // Constraints

            //- Apply constraints to equation
            template<class Type>
            void constrain(fvMatrix<Type>& eqn);

            //- Apply constraints to equation
            template<class Type>
            void constrain(fvBlockMatrix<Type>& eqn);

        // Correction

            //- Apply correction to field
            template<class Type>
            void correct(GeometricField<Type, fvPatchField, volMesh>& field);

            //- Correct fvOptions
            void correct();

            // auto-execute solver based on criterion
            void checkCorrect
            (
                const word& funcName,
                const word& fieldName
            );


        // Frame manipulations (encapsulates all MRF functions)

            //- Add the frame acceleration
            void addAcceleration
            (
                const volVectorField& U,
                volVectorField& ddtU
            ) const;

            //- Add the frame acceleration contribution to the momentum equation
            void addAcceleration
            (
                fvVectorMatrix& UEqn,
                bool rhs = false
            ) const;

            //- Add the frame acceleration contribution to the momentum equation
            void addAcceleration
            (
                const volScalarField& rho,
                fvVectorMatrix& UEqn,
                bool rhs = false
            ) const;

            //- Add the frame acceleration contribution to the momentum equation
            void addAcceleration
            (
                fvBlockMatrix<vector>& UEqn
            ) const;

            //- Add the frame acceleration contribution to the momentum equation
            void addAcceleration
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& UEqn
            ) const;

            //- Return the frame acceleration
            tmp<volVectorField> MRFDDt
            (
                const volVectorField& U
            );

            //- Return the frame acceleration
            tmp<volVectorField> MRFDDt
            (
                const volScalarField& rho,
                const volVectorField& U
            );


            //- Make the given absolute velocity relative within the MRF region
            void makeRelative(volVectorField& U) const;

            //- Make the given absolute flux relative within the MRF region
            void makeRelative(surfaceScalarField& phi) const;

            //- Return the given absolute flux relative within the MRF region
            tmp<surfaceScalarField> relative
            (
                const tmp<surfaceScalarField>& phi
            ) const;

            //- Return the given absolute boundary flux relative within
            //  the MRF region
            tmp<FieldField<fvsPatchField, scalar>> relative
            (
                const tmp<FieldField<fvsPatchField, scalar>>& tphi
            ) const;

            //- Return the given absolute patch flux relative within
            //  the MRF region
            tmp<Field<scalar>> relative
            (
                const tmp<Field<scalar>>& tphi,
                const label patchi
            ) const;

            //- Make the given absolute mass-flux relative within the MRF region
            void makeRelative
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Make the given relative velocity absolute within the MRF region
            void makeAbsolute(volVectorField& U) const;

            //- Make the given relative flux absolute within the MRF region
            void makeAbsolute(surfaceScalarField& phi) const;

            //- Return the given relative flux absolute within the MRF region
            tmp<surfaceScalarField> absolute
            (
                const tmp<surfaceScalarField>& phi
            ) const;

            //- Make the given relative mass-flux absolute within the MRF region
            void makeAbsolute
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Correct the boundary flux for the rotation of the MRF region
            void correctBoundaryFlux
            (
                const volVectorField& U,
                surfaceScalarField& phi
            ) const;

            tmp<surfaceScalarField> zeroFilter
            (
                const tmp<surfaceScalarField>& tphi
            ) const;

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

            //- Write data to Ostream
            virtual bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const optionList& options
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/fvOptions/fvOptionListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
