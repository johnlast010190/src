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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Class
    Foam::mappedFromFileFvPatchField

Description
    Foam::mappedFromFileFvPatchField

    A fixed value boundary mapped from an input file using a 0th order
    proximity method. The surface mesh is firstly read and its normals
    are constructed which allows distinction between front and back of
    simple baffles even if they are in the same patch.
    Although all file formats supported by the triSurface interfaces can
    be read in to create the surface mesh, only .nas and .ntl files are
    currently supported for the field values. Should be enhanced in
    the future.
    Not yet implemented: If a region or field is not present in the
    file, a warning will be printed and the provided value field
    will be used.

    Example of the boundary condition specification:
    @verbatim
    hotWall
    {
        type        mappedFromFile;
        file        "temperature.nas"; //relative to the case directory
        field       front_temp;
        convertToMillimeters true;
        value uniform 300; //required entry
    }
    @endverbatim


SourceFiles
    mappedFromFileFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFromFileFvPatchField_H
#define mappedFromFileFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/mappedFromFile/fileData.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mappedFromFileFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFromFileFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- data to map
        fileData<Type> fileData_;

        bool mappedData_;

    // Private member functions

        //- map get mapped fields from mapper
        tmp<Field<Type>> mapField
        (
            const pointField&,
            const vectorField&
        );

protected:

public:

    //- Runtime type information
    TypeName("mappedFromFile");


    // Constructors

        //- Construct from patch and internal field
        mappedFromFileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct by mapping given mappedFromFileFvPatchField
        //  onto a new patch
        mappedFromFileFvPatchField
        (
            const mappedFromFileFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct from patch, internal field and dictionary
        mappedFromFileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct as copy
        mappedFromFileFvPatchField
        (
            const mappedFromFileFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFromFileFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedFromFileFvPatchField
        (
            const mappedFromFileFvPatchField<Type>&,
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
                new mappedFromFileFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access


        // Mapping functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// Explicit specialisation declarations

template<>
tmp<Field<scalar>> mappedFromFileFvPatchField<scalar>::mapField
(
    const pointField&,
    const vectorField&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
