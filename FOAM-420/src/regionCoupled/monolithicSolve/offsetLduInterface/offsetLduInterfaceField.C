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
    (c) 2010-2015 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "offsetLduInterface/offsetLduInterfaceField.H"
#include "fields/Fields/Field/SubField.H"

#include "fields/fvPatchFields/constraint/processor/processorFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(offsetLduInterfaceField, 0);
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::offsetLduInterfaceField::offsetLduInterfaceField
(
    const lduInterfaceField& field,
    const lduInterface& slaveInterface,
    const lduInterface& offsetInterface,
    const label internalOffset,
    const label internalSize,
    const label foreignOffset,
    const label foreignSize
)
:
    lduInterfaceField(offsetInterface),
    slaveInterfaceField_(field),
    slaveInterface_(slaveInterface),
    internalOffset_(internalOffset),
    foreignOffset_(foreignOffset),
    internalSize_(internalSize),
    foreignSize_(foreignSize)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::offsetLduInterfaceField::~offsetLduInterfaceField()
{}


// * * * * * * * * * * * * * * Public members  * * * * * * * * * * * * * * * //

void Foam::offsetLduInterfaceField::initInterfaceMatrixUpdate
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction dirn,
    const Pstream::commsTypes commsType
) const
{
//    subResult_ = scalarField::subField(result, internalSize_, internalOffset_);
//timespec begin;
//clock_gettime(CLOCK_THREAD_CPUTIME_ID, &begin);
    scalarField::subField subResult(result, internalSize_, internalOffset_);
    slaveInterfaceField_.initInterfaceMatrixUpdate
    (
//        subResult_,
        const_cast<scalarField&>(static_cast<const scalarField&>(subResult)),
        add,
        scalarField::subField(psiInternal, foreignSize_, foreignOffset_),
        coeffs,
        dirn,
        commsType
    );
//timespec end;
//clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
//double dbgTime = (end.tv_nsec - begin.tv_nsec)*1e-9 + (end.tv_sec-begin.tv_sec);
//Info<< "offsetLduInterfaceField: " << dbgTime << endl;
//    SubList<scalar>(result, internalSize_, internalOffset_).assign(subResult_);
}

void Foam::offsetLduInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction dirn,
    const Pstream::commsTypes commsType
) const
{
    scalarField::subField subResult(result, internalSize_, internalOffset_);
//    subResult_ = scalarField::subField(result, internalSize_, internalOffset_);
    slaveInterfaceField_.updateInterfaceMatrix
    (
//        subResult_,
        const_cast<scalarField&>(static_cast<const scalarField&>(subResult)),
        add,
        scalarField::subField(psiInternal, foreignSize_, foreignOffset_),
        coeffs,
        dirn,
        commsType
    );
//    SubList<scalar>(result, internalSize_, internalOffset_).assign(subResult_);
//    subResult_.clear();
}


// ************************************************************************* //
