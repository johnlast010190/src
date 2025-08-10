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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::regIOobject::ownedByRegistry() const
{
    return ownedByRegistry_;
}


inline void Foam::regIOobject::store()
{
    if (!registered_)
    {
        FatalErrorInFunction
            << "Attempting to store an unregistered object: " << name() << nl
            << abort(FatalError);
    }
    ownedByRegistry_ = true;
}


template<class Type>
inline Type& Foam::regIOobject::store(Type* tPtr)
{
    if (!tPtr)
    {
        FatalErrorInFunction
            << "object deallocated"
            << abort(FatalError);
    }

    if (!tPtr->regIOobject::registered_)
    {
        FatalErrorInFunction
            << "Attempting to store an unregistered object." << nl
            << abort(FatalError);
    }
    tPtr->regIOobject::ownedByRegistry_ = true;

    return *tPtr;
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>& atPtr)
{
    Type* tPtr = atPtr.ptr();

    if (!tPtr)
    {
        FatalErrorInFunction
                << "object deallocated"
                << abort(FatalError);
    }

    if (!tPtr->regIOobject::registered_)
    {
        FatalErrorInFunction
            << "Attempting to store an unregistered object." << nl
            << abort(FatalError);
    }
    tPtr->regIOobject::ownedByRegistry_ = true;

    return *tPtr;
}


inline void Foam::regIOobject::release()
{
    ownedByRegistry_ = false;
}


inline Foam::label Foam::regIOobject::eventNo() const
{
    return eventNo_;
}

inline Foam::label& Foam::regIOobject::eventNo()
{
    return eventNo_;
}


inline const Foam::labelList& Foam::regIOobject::watchIndices() const
{
    return watchIndices_;
}


inline Foam::labelList& Foam::regIOobject::watchIndices()
{
    return watchIndices_;
}


// ************************************************************************* //
