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
    (c) 2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::patchMeanFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition sets its uniform fixed value from the mean of the
    same quantaty on another patch plus a specified off-set

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        patches      | source of uniform value | yes         | null
        offset       | offset from patch mean  | no          | zero
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            patchMean;
        patches        ( flowOut1 );

        offset constant 10;

        value uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    patchMeanFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef patchMeanFvPatchField_H
#define patchMeanFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class patchMeanFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class patchMeanFvPatchField
:
    public fixedValueFvPatchField<Type>
{

protected:

    // Protected data

        //- List of source patches to average
        List<wordRe> sourcePatchRegExp_;

        //- Corresponding patchIDs
        labelHashSet patchIDs_;

        //- Offset
        autoPtr<Function1<Type>> offset_;

public:

    //- Runtime type information
    TypeName("patchMean");


    // Constructors

        //- Construct from patch and internal field
        patchMeanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        patchMeanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchMeanFvPatchField
        //  onto a new patch
        patchMeanFvPatchField
        (
            const patchMeanFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        patchMeanFvPatchField
        (
            const patchMeanFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new patchMeanFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        patchMeanFvPatchField
        (
            const patchMeanFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new patchMeanFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Is defined in frame supported?
            virtual bool definedInFrameSupported() const
            {
                return true;
            }

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
