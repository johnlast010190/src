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
    (c) 2019-2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "algorithms/DelayedFunctionCall/DelayedFunctionCallImpl.H"


// * * * * * * * * * * * * * Static functions  * * * * * * * * * * * * * * * //

template <class ObjType, typename... ArgsType>
void Foam::serialThreads::run
(
    const label nActiveThreads,
    ObjType& obj,
    void (ObjType::*memFuncPtr)(ArgsType...),
    ArgsType&&... args
)
{
    if (!numThreads_)
    {
        FatalErrorInFunction
            << "run() called but a thread pool has not been created."
            << abort(FatalError);
    }
    if (nActiveThreads > numThreads_)
    {
        FatalErrorInFunction
            << "More threads requested than exist."
            << abort(FatalError);
    }

    functionCall_.set
    (
        new DelayedFunctionCallImpl
        <
            void,
            ObjType,
            void (ObjType::*)(ArgsType...),
            ArgsType&&...
        >
        (
            obj,
            memFuncPtr,
            std::forward<ArgsType>(args)...
        )
    );

    activeThread_ = -1;
    numActiveThreads_ = nActiveThreads;
    activateNext();
    // Master thread waits here while threads run
    wait();
    numActiveThreads_ = 0;
    functionCall_.clear();
}

// ************************************************************************* //
