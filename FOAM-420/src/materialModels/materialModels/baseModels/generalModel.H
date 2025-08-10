
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
    (c) 2022 Esi Ltd.

Class
    Foam::generalModel

\*---------------------------------------------------------------------------*/

#ifndef generalModel_H
#define generalModel_H

#include "fields/volFields/volFields.H"
#include "baseModels.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "materialModels/materialMacros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class generalModel Declaration
\*---------------------------------------------------------------------------*/

class generalModel
:
    public baseModels<scalar,fvPatchField,volMesh>
{
    // Private member data

        //- Changable model dimensions
        dimensionSet generalDims_;

        //- Changable function type
        word funcType_;

        //- Pointer to boundary nested class
        Boundary* boundary_;


public:

    TypeName("general");


    // Constructors

        //- Null constructor
        generalModel()
        {
            boundary_ = new Boundary(*this);
        };

        //- Destructor
        virtual ~generalModel()
        {
            delete boundary_;
        };

        //- Boundary class to handle boundaryField()[patchi] access
        class Boundary
        :
             public baseModels<scalar,fvPatchField,volMesh>::Boundary
        {
            public:
            Boundary(const generalModel& base)
            :
                generalModel_(base)
            {};
            const generalModel& generalModel_;
            virtual ~Boundary(){};
            virtual tmp<Field<scalar>> operator[](const label patchi) const
            {
                return generalModel_.generalPatch(patchi);
            }
        };


    // Public member functions

        //- Pure virtual function to compute geometric field
        virtual tmp<volScalarField> generalGeometric() const = 0;

        //- Pure virtual function to compute internal field
        virtual tmp<scalarField> generalInternal() const = 0;

        //- Pure virtual function to compute patch field
        virtual tmp<scalarField> generalPatch
        (
            label const patchi
        ) const = 0;

        //- Pure virtual function to compute cell values
        virtual scalar generalCell(label const celli) const = 0;

        //- Primitive field access/computation
        virtual tmp<scalarField> primitiveField() const override
        {
            return generalInternal();
        }

        //- Access to the boundary nested class for patch calculations
        virtual baseModels<scalar,fvPatchField,volMesh>::Boundary&
        boundaryField() const override
        {
            return *boundary_;
        }

        //- Model dimensions
        virtual const dimensionSet& dimensions() const override
        {
            return generalDims_;
        }

        //- Non-constant reference to allow dim overwrite
        virtual dimensionSet& dimensions()
        {
            return generalDims_;
        }


    // Operators

        //- Single cell access/computation
        virtual scalar operator[](const label celli) const override
        {
            return generalCell(celli);
        }

        //- Return geometric filed
        virtual tmp<GeometricField<scalar,fvPatchField,volMesh>>
        operator()() const override
        {
            return generalGeometric();
        }

    // Updates and additional functions

        //- Updating dependent model pointers
        virtual void updateTable(const word&) override = 0;

        //- Returns list of dependent model pointers
        virtual const depList& dep() const override = 0;

        //- Return function type
        virtual const word funcType() const override
        {
            return funcType_;
        }

        //- Type of the model
        virtual const word modelType() const override
        {
            return "scalar";
        }

        //- Pure virtual read function
        virtual bool read() override = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //