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
    (c) 2016-2017 Wikki Ltd
    (c) 2016, 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/constraint/processor/processorFaPatchField.H"
#include "faMesh/faPatches/constraint/processor/processorFaPatch.H"
#include "include/demandDrivenData.H"
#include "fields/Fields/transformField/transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{}


template<class Type>
processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    coupledFaPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)

{}


// Construct by mapping given processorFaPatchField<Type>
template<class Type>
processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)

{
    if (!isA<processorFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


template<class Type>
processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    coupledFaPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)

{
    if (!isA<processorFaPatch>(p))
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
}


template<class Type>
processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    procPatch_(refCast<const processorFaPatch>(ptf.patch())),
    sendBuf_(ptf.sendBuf_.xfer()),
    receiveBuf_(ptf.receiveBuf_.xfer()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(ptf.scalarSendBuf_.xfer()),
    scalarReceiveBuf_(ptf.scalarReceiveBuf_.xfer())
{
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}

template<class Type>
processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFaPatch>(ptf.patch())),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}



// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
processorFaPatchField<Type>::~processorFaPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> processorFaPatchField<Type>::patchNeighbourField() const
{
    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name()
            << " outstanding request."
            << abort(FatalError);
    }
    return *this;
}


template<class Type>
void processorFaPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        this->patchInternalField(sendBuf_);

        if (commsType == Pstream::commsTypes::nonBlocking && !Pstream::floatTransfer)
        {
            // Fast path. Receive into *this
            this->setSize(sendBuf_.size());
            outstandingRecvRequest_ = UPstream::nRequests();
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                procPatch_.neighbProcNo(),
                reinterpret_cast<char*>(this->begin()),
                this->byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );

            outstandingSendRequest_ = UPstream::nRequests();
            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                procPatch_.neighbProcNo(),
                reinterpret_cast<const char*>(sendBuf_.begin()),
                this->byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }
        else
        {
            procPatch_.compressedSend(commsType, sendBuf_);
        }
    }
}


template<class Type>
void processorFaPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if (commsType == Pstream::commsTypes::nonBlocking && !Pstream::floatTransfer)
        {
            // Fast path. Received into *this
            if
            (
                outstandingRecvRequest_ >= 0
             && outstandingRecvRequest_ < Pstream::nRequests()
            )
            {
                UPstream::waitRequest(outstandingRecvRequest_);
            }
            outstandingSendRequest_ = -1;
            outstandingRecvRequest_ = -1;
        }
        else
        {
            procPatch_.compressedReceive<Type>(commsType, *this);
        }

        if (doTransform())
        {
            Foam::transform(*this, procPatch_.forwardT(), *this);
        }
    }
}


template<class Type>
tmp<Field<Type>> processorFaPatchField<Type>::snGrad
(
        const scalarField& deltaCoeffs

) const
{
    return deltaCoeffs*(*this - this->patchInternalField());
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorFaPatchField<Type>::snGrad() const
{
    return this->patch().deltaCoeffs()*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
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

    const_cast<processorFaPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField&,
    const scalarField& coeffs,
    const direction cmpt,
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

        // Transform according to the transformation tensor
        transformCoupleField(scalarReceiveBuf_, cmpt);

        // Multiply the field by coefficients and add into the result
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

        // Transform according to the transformation tensor
        transformCoupleField(pnf, cmpt);

        // Multiply the field by coefficients and add into the result
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

    const_cast<processorFaPatchField<Type>&>(*this).updatedMatrix() = true;
}

template<class Type>
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, sendBuf_);

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


        receiveBuf_.setSize(sendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        IPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(receiveBuf_.begin()),
            receiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        OPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(sendBuf_.begin()),
            sendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, sendBuf_);
    }

    const_cast<processorFaPatchField<Type>&>(*this).updatedMatrix() = false;
}

template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>&,
    const scalarField& coeffs,
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

        // Consume straight from receiveBuf_

        // Transform according to the transformation tensor
        processorLduInterfaceField::transformCoupleField(receiveBuf_);

        // Multiply the field by coefficients and add into the result
//        forAll(edgeFaces, elemI)
//        {
//            result[edgeFaces[elemI]] -= coeffs[elemI]*receiveBuf_[elemI];
//        }
        if (!add)
        {
            forAll(edgeFaces, elemI)
            {
                result[edgeFaces[elemI]] += coeffs[elemI]*receiveBuf_[elemI];
            }
        }
        else
        {
            forAll(edgeFaces, elemI)
            {
                result[edgeFaces[elemI]] -= coeffs[elemI]*receiveBuf_[elemI];
            }
        }
    }
    else
    {
        Field<Type> pnf
        (
            procPatch_.compressedReceive<Type>(commsType, this->size())()
        );

        // Transform according to the transformation tensor
        processorLduInterfaceField::transformCoupleField(pnf);

        // Multiply the field by coefficients and add into the result
//        forAll(edgeFaces, elemI)
//        {
//            result[edgeFaces[elemI]] -= coeffs[elemI]*pnf[elemI];
//        }
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

    const_cast<processorFaPatchField<Type>&>(*this).updatedMatrix() = true;

}

template<class Type>
bool Foam::processorFaPatchField<Type>::ready() const
{
    if
    (
        outstandingSendRequest_ >= 0
     && outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingSendRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingSendRequest_ = -1;

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingRecvRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingRecvRequest_ = -1;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
