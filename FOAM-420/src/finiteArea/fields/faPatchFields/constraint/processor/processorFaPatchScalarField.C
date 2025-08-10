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
    (c) 2018 Esi Ltd.
    (c) 2016-2017 Wikki Ltd

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/constraint/processor/processorFaPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::processorFaPatchField<Foam::scalar>::transformCoupleField
(
    scalarField& f,
    const direction cmpt
) const
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::processorFaPatchField<Foam::scalar>::initInterfaceMatrixUpdate
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
//    procPatch_.compressedSend
//    (
//        commsType,
//        patch().patchInternalField(psiInternal)()
//    );
   this->patch().patchInternalField(psiInternal, scalarSendBuf_);

   if (commsType == Pstream::commsTypes::nonBlocking && !Pstream::floatTransfer)
   {
       // Fast path.
       if (debug && !this->ready())
       {
           FatalErrorInFunction
               << "On patch " << procPatch_.name()
               << " outstanding request."
               << abort(FatalError);
       }


       scalarReceiveBuf_.setSize(scalarSendBuf_.size());
       outstandingRecvRequest_ = UPstream::nRequests();
       UIPstream::read
       (
           Pstream::commsTypes::nonBlocking,
           procPatch_.neighbProcNo(),
           reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
           scalarReceiveBuf_.byteSize(),
           procPatch_.tag(),
           procPatch_.comm()
       );

       outstandingSendRequest_ = UPstream::nRequests();
       UOPstream::write
       (
           Pstream::commsTypes::nonBlocking,
           procPatch_.neighbProcNo(),
           reinterpret_cast<const char*>(scalarSendBuf_.begin()),
           scalarSendBuf_.byteSize(),
           procPatch_.tag(),
           procPatch_.comm()
       );
   }
   else
   {
       procPatch_.compressedSend(commsType, scalarSendBuf_);
   }

   const_cast<processorFaPatchField<scalar>&>(*this).updatedMatrix() = false;

}


template<>
void processorFaPatchField<scalar>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField&,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& edgeFaces = this->patch().edgeFaces();

    if (commsType == Pstream::commsTypes::nonBlocking && !Pstream::floatTransfer)
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }
        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;


        // Consume straight from scalarReceiveBuf_
        if (!add)
        {
            forAll(edgeFaces, elemI)
            {
                result[edgeFaces[elemI]] += coeffs[elemI]*scalarReceiveBuf_[elemI];
            }
        }
        else
        {
            forAll(edgeFaces, elemI)
            {
                result[edgeFaces[elemI]] -= coeffs[elemI]*scalarReceiveBuf_[elemI];
            }
        }
    }
    else
    {
        scalarField pnf
        (
            procPatch_.compressedReceive<scalar>(commsType, this->size())()
        );

        if (!add)
        {
            forAll(edgeFaces, elemI)
            {
                result[edgeFaces[elemI]] += coeffs[elemI]*pnf[elemI];
            }
        }
        else
        {
            forAll(edgeFaces, elemI)
            {
                result[edgeFaces[elemI]] -= coeffs[elemI]*pnf[elemI];
            }
        }
    }

    const_cast<processorFaPatchField<scalar>&>(*this).updatedMatrix() = true;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
