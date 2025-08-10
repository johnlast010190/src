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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/lduAddressing/lduInterface/processorLduInterface.H"
#include "db/IOstreams/Pstreams/IPstream.H"
#include "db/IOstreams/Pstreams/OPstream.H"

// * * * * * * * * * * * * * * * Member Functions * * *  * * * * * * * * * * //

template<class Type>
void Foam::processorLduInterface::send
(
    const Pstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    label nBytes = f.byteSize();

    if
    (
        commsType == Pstream::commsTypes::blocking
     || commsType == Pstream::commsTypes::scheduled
    )
    {
        OPstream::write
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<const char*>(f.begin()),
            nBytes,
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::commsTypes::nonBlocking)
    {
        resizeBuf(receiveBuf_, nBytes);

        IPstream::read
        (
            commsType,
            neighbProcNo(),
            receiveBuf_.begin(),
            nBytes,
            tag(),
            comm()
        );

        resizeBuf(sendBuf_, nBytes);
        memcpy(sendBuf_.begin(), f.begin(), nBytes);

        OPstream::write
        (
            commsType,
            neighbProcNo(),
            sendBuf_.begin(),
            nBytes,
            tag(),
            comm()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << exit(FatalError);
    }
}


template<class Type>
void Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if
    (
        commsType == Pstream::commsTypes::blocking
     || commsType == Pstream::commsTypes::scheduled
    )
    {
        IPstream::read
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<char*>(f.begin()),
            f.byteSize(),
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::commsTypes::nonBlocking)
    {
        memcpy
        (
            reinterpret_cast<char*>(f.begin()),
            receiveBuf_.begin(),
            f.byteSize()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << exit(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    const label size
) const
{
    tmp<Field<Type>> tf(new Field<Type>(size));
    receive(commsType, tf.ref());
    return tf;
}


template<class Type>
void Foam::processorLduInterface::compressedSend
(
    const Pstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        const scalar *sArray = reinterpret_cast<const scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];
        resizeBuf(sendBuf_, nBytes);
        float *fArray = reinterpret_cast<float*>(sendBuf_.begin());

        for (label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        // Type might be a type with greater alignment requirements than float, so use a memcpy.
        memcpy(sendBuf_.begin() + nm1 * sizeof(float), &f.last(), sizeof(Type));

        if
        (
            commsType == Pstream::commsTypes::blocking
         || commsType == Pstream::commsTypes::scheduled
        )
        {
            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType == Pstream::commsTypes::nonBlocking)
        {
            resizeBuf(receiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                receiveBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );

            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported communications type " << int(commsType)
                << exit(FatalError);
        }
    }
    else
    {
        this->send(commsType, f);
    }
}

template<class Type>
void Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        if
        (
            commsType == Pstream::commsTypes::blocking
         || commsType == Pstream::commsTypes::scheduled
        )
        {
            resizeBuf(receiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                receiveBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType != Pstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Unsupported communications type " << int(commsType)
                << exit(FatalError);
        }

        memcpy(&f.last(), receiveBuf_.begin() + nm1 * sizeof(float), sizeof(Type));

        const float *fArray =
            reinterpret_cast<const float*>(receiveBuf_.begin());
        scalar *sArray = reinterpret_cast<scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];

        for (label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }
    }
    else
    {
        this->receive<Type>(commsType, f);
    }
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    const label size
) const
{
    tmp<Field<Type>> tf(new Field<Type>(size));
    compressedReceive(commsType, tf.ref());
    return tf;
}


// ************************************************************************* //
