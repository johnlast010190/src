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
    Foam::referenceFrameFvPatch

Description
    Class that helps manipulation of boundary conditions when they are bounded
    to a referenceFrame.
\*---------------------------------------------------------------------------*/

#ifndef referenceFrameFvPatch_H
#define referenceFrameFvPatch_H

#include "fvMesh/fvMesh.H"
#include "primitives/bools/Switch/Switch.H"
#include "referenceFrames/coordinateFrame.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class referenceFrameFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class referenceFrameFvPatch
{
    // Private memeber funcitons

        //- Helper function to load name of U field when possible
        word loadUName(const DimensionedField<Type, volMesh>& iF) const;

        //- Helper function to load object registry
        const objectRegistry& loadDb
        (
            const DimensionedField<Type, volMesh>& iF
        ) const;


protected:

    // Protected member data

        //- Assignable coordinate frame
        mutable coordinateFrame* coorFramePtr_;

        //- Reference to the patch
        const fvPatch& patch_;

        //- Sotrage for the input data
        Field<Type> inputValue_;

        //- Is inlet flux on?
        Switch inletFlux_;

        //- Name of reference frame
        word frameName_;

        //- Velocity field name
        word refPatchUName_;

        //- Reference to object registry
        const objectRegistry& obr_;

        //- Is the patch data defined in frame?
        Switch definedInFrame_;


    // Protected member functions

        //- Setting the coordinate frame pointer
        void setCoorFramePtr() const;


public:

    // Constructors

        //- Construct from patch
        explicit referenceFrameFvPatch
        (
            const fvPatch& patch,
            const DimensionedField<Type, volMesh>& iF =
                NullObjectRef<DimensionedField<Type, volMesh>>()
        );

        //- Construct from data
        referenceFrameFvPatch
        (
            const fvPatch& patch,
            coordinateFrame* coorFramePtr,
            const Field<Type>& inputValue,
            Switch inletFlux,
            const word& frameName,
            const DimensionedField<Type, volMesh>& iF =
                NullObjectRef<DimensionedField<Type, volMesh>>()
        );

        //- Construct from data and mapper
        referenceFrameFvPatch
        (
            const fvPatch& patch,
            coordinateFrame* coorFramePtr,
            const Field<Type>& inputValue,
            Switch inletFlux,
            const word& frameName,
            const fvPatchFieldMapper& mapper,
            const DimensionedField<Type, volMesh>& iF =
                NullObjectRef<DimensionedField<Type, volMesh>>()
        );

        //- Construct from dict and patch
        referenceFrameFvPatch
        (
            const dictionary& dict,
            const fvPatch& patch,
            const DimensionedField<Type, volMesh>& iF =
                NullObjectRef<DimensionedField<Type, volMesh>>()
        );


    //- Destructor
    virtual ~referenceFrameFvPatch() = default;


    // Member functions

        //- Return pointer to coordinate frame
        const coordinateFrame* coorFramePtr() const
        {
            if (!coorFramePtr_)
            {
                setCoorFramePtr();
            }
            return coorFramePtr_;
        }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap
            (
                const referenceFrameFvPatch<Type>& ptf,
                const labelList& addr
            );

            //- Map (and resize as needed) from self given a mapping object
            void autoMapGIB(const gibFvPatchFieldMapper& mapper);


        // Evaluation functions

            //- Function adds reference frame velocity
            // It optionaly updates stored field
            void addFrameVelocity
            (
                Field<Type>& Up,
                bool setInValue = true,
                bool setToGlobal = true
            );

            //- Function adds reference frame gradients
            void addFrameGradient(Field<Type>& gradient) const;

            //- Local to global transformation
            void makeVectorGlobal(Field<Type>& U) const;

            //- Local to global transformation
            void makeVectorGlobal(Type& U) const;

            //- Making fluxes relative
            void makeRelative(scalarField& phiP);

            //- Making the patch velocity relative
            void makeRelative(vectorField& Up);

            //- Attach coordinate frame to patch
            void updateCordinateFrameRegistry() const;

            //- Returns frame velocity
            tmp<Field<Type>> getFrameVelocity() const;

            //- There are cases (like totalPressure boundary conditions) that we
            // don't want to use the templated version of the getFrameVelocity
            tmp<vectorField> getFrameVelocityVectorField();

            //- Name of the reference frame
            const word& frameName() const
            {
                return frameName_;
            }

            //- Is inlet flux on?
            Switch inletFlux() const
            {
                return inletFlux_;
            }

            //- Reference to input value
            const Field<Type>& inputValue() const
            {
                return inputValue_;
            }

            //- Re-setting input value
            void setInputValue(const Field<Type>& ft)
            {
                inputValue_ = ft;
            }

            //- Is defined in frame supported?
            virtual bool definedInFrameSupported() const
            {
                if (definedInFrame_)
                {
                    FatalErrorInFunction
                        << "The boundary condition for patch: "
                        << "\"" << patch_.name() << "\""
                        << " field: \"" << refPatchUName_ << "\""
                        << " doesn't support definition "
                        << "within the frame." << exit(FatalError);
                }
                return false;
            }

            //- Is defined in frame?
            virtual bool isDefinedInFrame() const
            {
                return definedInFrame_;
            }


    //- Read reference frame info
    void read(const dictionary& dict);

    //- Write settings
    void write(Ostream&) const;
};


template<>
void referenceFrameFvPatch<vector>::addFrameVelocity
(
    vectorField& Up,
    bool setInValue,
    bool setToGlobal
);


template<>
void referenceFrameFvPatch<vector>::addFrameGradient
(
    vectorField& gradient
) const;


template<>
void referenceFrameFvPatch<vector>::makeVectorGlobal(vectorField& U) const;


template<>
void referenceFrameFvPatch<vector>::makeVectorGlobal(vector& U) const;


template<>
tmp<vectorField>
referenceFrameFvPatch<vector>::getFrameVelocity() const;


template<>
inline tmp<scalarField>
referenceFrameFvPatch<scalar>::getFrameVelocity() const
{
    return tmp<scalarField>(new scalarField(patch_.size(), 0.));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
