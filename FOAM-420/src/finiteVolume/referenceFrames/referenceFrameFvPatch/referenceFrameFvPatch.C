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
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatchTemplates.C"
#include "finiteVolume/fvc/fvcMeshPhi.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::referenceFrameFvPatch<Foam::vector>::addFrameVelocity
(
    vectorField& Up,
    bool setInValue,
    bool setToGlobal
)
{
    if (setInValue)
    {
        setInputValue(Up);
    }

    // Transforming velocity from the local frame velocity to global
    if (setToGlobal)
    {
        makeVectorGlobal(Up);
    }

    if (coorFramePtr())
    {
        coorFramePtr()->calculateBoundaryVelocity
        (
            Up,
            patch_.index(),
            obr_,
            refPatchUName_,
            inletFlux_
        );
    }
}


template<>
void Foam::referenceFrameFvPatch<Foam::vector>::addFrameGradient
(
    vectorField& gradient
) const
{
    // Transforming velocity from the local frame velocity to global
    makeVectorGlobal(gradient);

    // So far the gradient is just rotated in the frame but other forms
    // of gradients could be implemented comming from the frame.
}


template<>
void Foam::referenceFrameFvPatch<Foam::vector>::makeVectorGlobal
(
    vectorField& U
) const
{
    if (coorFramePtr() && definedInFrameSupported() && isDefinedInFrame())
    {
        U = coorFramePtr()->coorSys().globalVector(U);
    }
}


template<>
void Foam::referenceFrameFvPatch<Foam::vector>::makeVectorGlobal
(
    vector& U
) const
{
    if (coorFramePtr() && definedInFrameSupported() && isDefinedInFrame())
    {
        U = coorFramePtr()->coorSys().globalVector(U);
    }
}


template<>
Foam::tmp<Foam::Field<Foam::vector>>
Foam::referenceFrameFvPatch<Foam::vector>::getFrameVelocity() const
{
    tmp<vectorField> Uframe(new vectorField(patch_.size(), Zero));
    if (coorFramePtr())
    {
        coorFramePtr()->calculateBoundaryVelocity
        (
            Uframe.ref(),
            patch_.index(),
            obr_,
            refPatchUName_,
            true
        );
        return Uframe;
    }
    return Uframe;
}


// * * * * * * * * * * * * * * * Instantiations  * * * * * * * * * * * * * * //

namespace Foam
{
    INSTANTIATE_FOR_ALL_FIELD_TYPES(referenceFrameFvPatch);
    INSTANTIATE_FOR_ALL_VECTOR_TENSOR_N_TYPES(referenceFrameFvPatch);
}

// ************************************************************************* //
