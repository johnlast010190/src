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
    (c) 2020-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIFvPatchVectorNFields_H
#define cyclicAMIFvPatchVectorNFields_H

#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"
#include "VectorN/finiteVolume/fields/fvPatchFields/fvPatchVectorNFields.H"

#include "fields/fvPatchFields/constraint/cyclicAMI/cyclicAMIFvPatchFields.H"
#include "VectorN/finiteVolume/fields/fvPatchFields/cyclicAMIFvPatchVectorNFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeVectorNCyclic(Type)                                                \
template <>                                                                    \
tmp<Field<Type>> cyclicAMIFvPatchField<Type>::patchNeighbourField              \
(                                                                              \
    const Pstream::commsTypes                                                  \
) const;                                                                       \
                                                                               \
template <>                                                                    \
tmp<Field<Type>> cyclicAMIFvPatchField<Type>::patchInternalField() const;     \
                                                                               \
template <>                                                                    \
void cyclicAMIFvPatchField<Type>::initInterfaceMatrixUpdate                    \
(                                                                              \
    const Field<Type>& psiInternal,                                            \
    Field<Type>&,                                                              \
    const BlockLduMatrix<Type>&,                                               \
    const CoeffField<Type>&,                                                   \
    const Pstream::commsTypes commsType                                        \
) const;                                                                       \
                                                                               \
template <>                                                                    \
void cyclicAMIFvPatchField<Type>::updateInterfaceMatrix                        \
(                                                                              \
    const Field<Type>& psiInternal,                                            \
    Field<Type>& result,                                                       \
    const BlockLduMatrix<Type>&,                                               \
    const CoeffField<Type>& coeffs,                                            \
    const Pstream::commsTypes commsType                                        \
) const;\
                                                                               \
template <>                                                                    \
void cyclicAMIFvPatchField<Type>::addToInternalCoupledField                    \
(                                                                              \
    Field<Type>&,                                                              \
    const Field<Type>&,                                                        \
    const CoeffField<Type>&                                                    \
) const;\
\
template <>                                                                    \
void cyclicAMIFvPatchField<Type>::addToInternalCoupledCoveredField             \
(                                                                              \
    Field<Type>&,                                                              \
    const Field<Type>&,                                                        \
    const CoeffField<Type>&                                                    \
) const;


#define doMakePatchTypeFieldDef(type, Type, args...)                           \
    makeVectorNCyclic(type)

forAllVectorNTypes(doMakePatchTypeFieldDef)

#undef doMakePatchTypeFieldDef

#undef makeVectorNCyclic

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
