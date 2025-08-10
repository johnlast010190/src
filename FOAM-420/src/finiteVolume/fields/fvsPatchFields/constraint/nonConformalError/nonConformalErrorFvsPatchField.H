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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::nonConformalErrorFvsPatchField

Description
    Constraint fvsPatchField for nonConformalError couplings.

SourceFiles
    nonConformalErrorFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorFvsPatchField_H
#define nonConformalErrorFvsPatchField_H

#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalError/nonConformalErrorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nonConformalErrorFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalErrorFvsPatchField
:
    public fvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(nonConformalErrorFvPatch::typeName_());


    // Constructors

        //- Construct from patch and nonConformalError field
        nonConformalErrorFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, nonConformalError field and dictionary
        nonConformalErrorFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalErrorFvsPatchField onto a
        //  new patch
        nonConformalErrorFvsPatchField
        (
            const nonConformalErrorFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting nonConformalError field reference
        nonConformalErrorFvsPatchField
        (
            const nonConformalErrorFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting nonConformalError field reference
        nonConformalErrorFvsPatchField
        (
            const nonConformalErrorFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting nonConformalError field
        //  reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const override
        {
            return tmp<fvsPatchField<Type>>
            (
                new nonConformalErrorFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvsPatchFields/constraint/nonConformalError/nonConformalErrorFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
