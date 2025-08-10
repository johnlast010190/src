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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.


\*---------------------------------------------------------------------------*/

#include "faMesh/faPatches/constraint/wedge/wedgeFaPatch.H"
#include "fields/faPatchFields/constraint/wedge/wedgeFaPatchField.H"
#include "fields/Fields/transformField/transformField.H"
#include "primitives/transform/symmTransform.H"
#include "primitives/DiagTensor/diagTensor/diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    if (!isType<wedgeFaPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    this->evaluate();
}


template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return gradient at boundary
template<class Type>
tmp<Field<Type>> wedgeFaPatchField<Type>::snGrad() const
{
    const Field<Type> pif(this->patchInternalField());

    return
    (
        transform(refCast<const wedgeFaPatch>(this->patch()).faceT(), pif)
      - pif
    )*(0.5*this->patch().deltaCoeffs());
}


// Evaluate the patch field
template<class Type>
void wedgeFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    faPatchField<Type>::forceAssign
    (
        transform
        (
            refCast<const wedgeFaPatch>(this->patch()).edgeT(),
            this->patchInternalField()
        )
    );
}


// Return defining fields
template<class Type>
tmp<Field<Type>> wedgeFaPatchField<Type>::snGradTransformDiag() const
{
    const diagTensor diagT =
        0.5*diag(I - refCast<const wedgeFaPatch>(this->patch()).faceT());

    const vector diagV(diagT.xx(), diagT.yy(), diagT.zz());

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            this->size(),
            transformMask<Type>
            (
                pow
                (
                    diagV,
                    pTraits<typename powProduct<vector, pTraits<Type>::rank>
                    ::type>::zero
                )
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
