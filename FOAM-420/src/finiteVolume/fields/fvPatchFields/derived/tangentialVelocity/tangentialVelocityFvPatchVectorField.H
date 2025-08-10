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
    (c) 2010-2012 Esi Ltd.
    (c) 2006-2008 OpenCFD Ltd.

Class
    Foam::tangentialVelocityFvPatchVectorField

Description
    A boundary for providing tangential velocity on walls. The boundary supports
    two input modes: the user can either supply a "cross" vector or a direction
    vector.

    When a cross vector is supplied, the velocity is calculated from:

    Direction = (Cross ^ normal)/(mag(Cross ^ normal)

    U = Direction* Umag - normal * ((normal & Direction) * Umag)

    When a Direction is supplied, the velocity is simply calculated from:

    U = Direction* Umag - normal * ((normal & Direction) * Umag)

    If both Direction and Cross vector are supplied, the Cross vector takes
    precedence.

    Example of the boundary condition specification:
    @verbatim
    roller
    {
        type            tangentialVelocity;
        magnitude       10; //velocity magnitude on properly aligned faces
        crossVector     (0 1 0); //Cross direction
        direction       (1 0 0); //not used if crossVector exists
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

    Notes:

    -Faces not parallel to the direction vector will have a velcoity
    magnitude smaller than the specified magnitude.

    -The crossVector entry is required if used in conjunction with the
    tangentialBeltVelocity boundary



SourceFiles
    tangentialVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tangentialVelocityFvPatchVectorField_H
#define tangentialVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class tangentialVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tangentialVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- magnitude of the wall velocity
        scalar magVel_;

        //- the vector defining third direction.
        //-tangential velocity will be this cross normal vector
        autoPtr<vector> crossVectorPtr_;

        //- Vector defining velocity direction before
        //- non-tangential components are removed
        autoPtr<vector> directionPtr_;


public:

   //- Runtime type information
   TypeName("tangentialVelocity");


   // Constructors

        //- Construct from patch and internal field
        tangentialVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tangentialVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  tangentialVelocityFvPatchVectorField
        //  onto a new patch
        tangentialVelocityFvPatchVectorField
        (
            const tangentialVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tangentialVelocityFvPatchVectorField
        (
            const tangentialVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tangentialVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tangentialVelocityFvPatchVectorField
        (
            const tangentialVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tangentialVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
        // Access

            //- Return the rotation orientation vector
            vector crossVector() const
            {
                if (crossVectorPtr_.operator->())
                {
                    return crossVectorPtr_();
                }

                FatalErrorInFunction
                    << "Call for undefined variable "
                    << "crossVectorPtr_ in patch "
                    << this->patch().name()
                    << exit(FatalError);
                return crossVectorPtr_();
            }


            //- Return the velocity magnitude vector
            scalar magVel() const
            {
                return magVel_;
            }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Is defined in frame supported?
        virtual bool definedInFrameSupported() const
        {
            return true;
        }

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
