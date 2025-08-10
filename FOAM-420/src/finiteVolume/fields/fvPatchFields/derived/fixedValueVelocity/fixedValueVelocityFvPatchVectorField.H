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
    (c) 2019 Esi Ltd.

Class
    Foam::fixedValueVelocityFvPatchVectorField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a fixed value velocity which is aware
    of the reference frame that is being assigned to.

Usage

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedValaueVelocity;
        referenceFrame  frameName;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::coordinateFrame

SourceFiles
    fixedValueVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueVelocityFvPatchVectorField_H
#define fixedValueVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class fixedValueVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedValueVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("fixedValueVelocity");


    // Constructors

        //- Construct from patch and internal field
        fixedValueVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueVelocityFvPatchVectorField
        //  onto a new patch
        fixedValueVelocityFvPatchVectorField
        (
            const fixedValueVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueVelocityFvPatchVectorField
        (
            const fixedValueVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new fixedValueVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedValueVelocityFvPatchVectorField
        (
            const fixedValueVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        virtual ~fixedValueVelocityFvPatchVectorField()
        {}

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
             const DimensionedField<vector, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchVectorField>
            (
                new fixedValueVelocityFvPatchVectorField(*this, iF)
            );
        }

    // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        ) override;

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchVectorField&,
            const labelList&
        ) override;

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs() override;

        //- Is defined in frame supported?
        virtual bool definedInFrameSupported() const override
        {
            return true;
        }

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

