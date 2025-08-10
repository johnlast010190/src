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
    (c) 2018 Esi Ltd.

Class
    Foam::interpolatedInletOutletFvPatchField


Description
    This boundary condition provides spatially 1D varying inletOutlet.
    Inflow profile is specified via Function1
    X-coordinate is supplied via patchDistanceFunction
     - Wall distance
     - Distance along a vector
     - X, Y, Z coordinate
     - Distance from point

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        profile      | profile entry           | yes         |
        distanceType | distance function       | yes         |
        xscale       | scaling for distance    | no          | 1
        yscale       | scaling for profile     | no          | 1
        xoffset      | distance additive const | no          | 0
        yoffset      | profile additive const  | no          | 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            interpolatedInletOutlet;
        profile         table
        (
            (0  1)
            (1  1.5)
            (5  2)
        );
        distanceType    wallDistance;
        xscale          1;
        yscale          1;
        xoffset         0;
        yoffset         0;

        value uniform 0;
    }
    \endverbatim

Note
    This class is derived directly from a mixed patch to allow inletOut type
    behaviour

See Also
    Foam::interpolationTable and Foam::mixedFvPatchField

SourceFiles
    interpolatedInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedInletOutletFvPatchField_H
#define interpolatedInletOutletFvPatchField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchField.H"
#include "cfdTools/general/fieldProfile/fieldProfile.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class interpolatedInletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolatedInletOutletFvPatchField
:
    public mixedFvPatchField<Type>
{
private:
    // Private data

        //- The distribution being used, including the bounding treatment
        autoPtr<fieldProfile<Type>> distribution_;

        //- flux name, "none" for fixedValue
        word phiName_;


public:

    //- Runtime type information
    TypeName("interpolatedInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        interpolatedInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interpolatedInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        interpolatedInletOutletFvPatchField
        (
            const interpolatedInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interpolatedInletOutletFvPatchField
        (
            const interpolatedInletOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new interpolatedInletOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interpolatedInletOutletFvPatchField
        (
            const interpolatedInletOutletFvPatchField<Type>&,
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
                new interpolatedInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate dY/dX for face centres
            tmp<Field<Type>> gradient () const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
