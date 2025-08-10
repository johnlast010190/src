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
    (c) 2010 Ivor Clifford

\*---------------------------------------------------------------------------*/

#include "VectorN/finiteVolume/fields/fvPatchFields/wedgeFvPatchVectorNFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/constraint/wedge/wedgeFvPatchField.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeVectorTensorNWedgeFuncDefs(Type)                    \
template<>                                                      \
tmp<Field<Type>> wedgeFvPatchField<Type>::snGrad() const       \
{                                                               \
    return tmp<Field<Type>>                                    \
    (                                                           \
        new Field<Type>(size(), pTraits<Type>::zero)            \
    );                                                          \
}                                                               \
                                                                \
template<>                                                      \
void wedgeFvPatchField<Type>::evaluate(                         \
    const Pstream::commsTypes commsType                         \
)                                                               \
{                                                               \
    if (!updated())                                             \
    {                                                           \
        updateCoeffs();                                         \
    }                                                           \
                                                                \
    wedgeFvPatchField<Type>::forceAssign(patchInternalField());  \
}                                                               \
                                                                \
template<>                                                      \
tmp<Field<Type>> wedgeFvPatchField<Type>::snGradTransformDiag()\
const                                                           \
{                                                               \
    return tmp<Field<Type>>                                    \
    (                                                           \
        new Field<Type>(this->size(), pTraits<Type>::zero)      \
    );                                                          \
}


#define doMakePatchTypeField(type, Type, args...)                           \
    makeVectorTensorNWedgeFuncDefs(type)                                    \
                                                                            \
    makeTemplatePatchTypeField(type, wedge);


forAllVectorNTypes(doMakePatchTypeField)

forAllTensorNTypes(doMakePatchTypeField)

forAllDiagTensorNTypes(doMakePatchTypeField)

forAllSphericalTensorNTypes(doMakePatchTypeField)


#undef doMakePatchTypeField

#undef makeVectorTensorNWedgeFuncDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
