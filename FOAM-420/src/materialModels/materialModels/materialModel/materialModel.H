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
    (c) 2021-2023 Esi Ltd.

Class
    Foam::materialModel

Description
    Base class for selecting material models.

SourceFiles
    materialModel.C
    materialModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef materialModel_H
#define materialModel_H

#include "materialModels/materialTables/materialTables.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class materialModel Declaration
\*---------------------------------------------------------------------------*/

class materialModel
:
    public regIOobject
{
    // Private member functions

        //- Return reference to the mesh
        const fvMesh& meshFromRegistry
        (
            const objectRegistry& obr
        ) const;


protected:

    // Protected data

        //- Materials object registry
        const objectRegistry& matObr_;

        //- Material models tables
        materialTables& materialTables_;

        //- Phase name
        const word phaseName_;

        //- Specie name
        const word specieName_;

        //- Model name (in object registry)
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Object registry reference
        const objectRegistry& obr_;

        //- Pointer to dictionary/sub-dictionary
        const dictionary* dict_;

        //- Dependance ptr list
        depList dep_;

        //- Pointer list to the models
        matScalarUTable sMod_;

        //- Pointer list to the models
        matVectorUTable vMod_;

        //- Pointer list to the models
        matTensorUTable tMod_;

        //- Pointer to the temperature offset by reference
        const refScalarField* T_;

        //- Pointer to the pressure offset by the reference
        const refScalarField* p_;


    // Protected member functions

        //- Look up an optional dependent model and return if it is constant
        template<class Type>
        bool isModelConst() const;


public:

    //- Runtime type information
    TypeName("materialModel");

    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        materialModel,
        dictionary,
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        ),
        (obr, dict, phaseName, specieName, name)
    );


    // Selectors

        //- Return reference to selected rho model
        static autoPtr<materialModel> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& modelType,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Constructors

        //- Construct from components
        materialModel
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name = word::null
        );


    //- Destructor
    virtual ~materialModel();


    // Public member functions

        //- Constructs phase filed name
        static word phasePropertyName
        (
            const word fieldName,
            const word& phaseName = word::null
        );

        //- Fill dependencies
        //  First check is if checkModel == modelName
        //  Than dependencies from modelName are added
        void fill
        (
            const word& checkModel,         // Check if the model name corresponds to this
            const word& modelName,          // Model for which dependencies are filled
            const wordList& modelNames,     // Required models for modelName
            const labelList& inds,          // Indexes to put dependent models in
            const boolList& modelsRequired, // Each model has a flag if it is required or not
            const label depNumber           // Slot number for overall dep_ list
        );

        //- Return first cell value from model
        scalar value
        (
            const word& modelName,
            const label modelIndex,
            const label celli = 0,
            const bool update = true,
            const bool readDict = false
        );

        //- Constructs phase filed name
        void updateDictPtr();

        const dictionary& materialsDict() const;

        //- Get field pointer with the nullptr check
        template<class Type>
        GeometricField<Type, fvPatchField, volMesh>* lookupPtr
        (
            const word& fieldSubName,
            const word& phaseName = word::null
        );

        //- Constructiong or returning pointer to reference field
        template<class Type>
        const referenceFields<Type>* constructOrReturnRefFieldPtr
        (
            const word& fieldName
        ) const;

        //- Mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Update table
        virtual void updateTable(const word&) = 0;

        //- Type independent pointer list
        virtual const depList& dep() const
        {
            return dep_;
        }

        //- Casting of scalar models
        virtual baseModels<scalar>* castScalarModel(const word&)
        {
            return nullptr;
        }

        //- Casting of vector models
        virtual baseModels<vector>* castVectorModel(const word&)
        {
            return nullptr;
        }

        //- Casting of tensor models
        virtual baseModels<tensor>* castTensorModel(const word&)
        {
            return nullptr;
        }

        //- Explicitly update scalar field
        virtual bool updateScalarField
        (
            const word& fieldName,
            const volScalarField& volField
        )
        {
            return false;
        };

        //- Explicitly update vector field
        virtual bool updateVectorField
        (
            const word& fieldName,
            const volVectorField& volField
        )
        {
            return false;
        };

        //- Explicitly update tensor field
        virtual bool updateTensorField
        (
            const word& fieldName,
            const volTensorField& volField
        )
        {
            return false;
        };

        //- Is it mixture model
        virtual bool mixture() const
        {
            return false;
        }

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return false;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const
        {
            return false;
        }

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return true;
        }

        //- Return true if Cp or Cv is uniform depending on energy type
        virtual bool isCpvConst() const
        {
            return false;
        }

        //- Write for data
        virtual bool writeData(Ostream& os) const override
        {
            NotImplemented;
        }

        //- Write
        virtual bool write(const bool valid) const override
        {
            NotImplemented;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "materialModelTemplates.C"
#endif

#endif


// ************************************************************************* //
