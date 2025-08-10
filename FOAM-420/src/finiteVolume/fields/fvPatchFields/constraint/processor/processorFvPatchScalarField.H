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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2011 Symscape
    (c) 2014 Symscape
    (c) 2019 Esi Ltd.

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.

\*---------------------------------------------------------------------------*/

#ifndef processorFvPatchScalarField_H
#define processorFvPatchScalarField_H

#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void processorFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    scalarField&,
    const bool add,
    const scalarField&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const;


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField&,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const;


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    const Field<scalar>& psiInternal,
    Field<scalar>& result,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>& coeffs,
    const Pstream::commsTypes commsType
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
