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
    (c) 2010-2011, Esi Ltd

Class
    Foam::resistiveVelocityFvPatchVectorField

Description
    Set the velocity normal gradient for a resistive outlet
    Normal velocity: Un = Uin
    Tangential velocity: Ut = Ct * Uit (Ct from pressure boundary)
    Where Ui is the patchInternal velocity field

    Example of the boundary condition specification:
    @verbatim
    outlet
    {
        type            resistiveVelocity;
        //p               p;
        //phi              phi;
        //value           uniform (0 0 0); // placeholder/result/restart
        //gradient        uniform (0 0 0); // placeholder/result/restart
    }
    @endverbatim

SourceFiles
    resistiveVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef resistiveVelocityFvPatchVectorFields_H
#define resistiveVelocityFvPatchVectorFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/directionMixed/directionMixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class resistiveVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class resistiveVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- pressure field name
        word pName_;

        //- flux field name
        word phiName_;

public:

    //- Runtime type information
    TypeName("resistiveVelocity");


    // Constructors

        //- Construct from patch and internal field
        resistiveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        resistiveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  resistiveVelocityFvPatchVectorField onto a new patch
        resistiveVelocityFvPatchVectorField
        (
            const resistiveVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        resistiveVelocityFvPatchVectorField
        (
            const resistiveVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new resistiveVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        resistiveVelocityFvPatchVectorField
        (
            const resistiveVelocityFvPatchVectorField&,
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
                new resistiveVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- use weak boundary type for tangential velocity
        //  set value coefficient to zero for outflow
        virtual tmp<Field<vector>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
