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
    (c) 2021 Esi Ltd.

Namespace
    parallelForAll

Description
    Tool for code parallelization with TBB

\*---------------------------------------------------------------------------*/

#ifndef ForEach_H
#define ForEach_H

#include "stdFoam.H"

#ifdef FOAM_USE_TBB
    #include <tbb/parallel_for.h>
    #include <tbb/parallel_for_each.h>
    #include <tbb/blocked_range.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    template<typename T, typename Func>
    void parallelForAll
    (
        const UList<T>& list,
        const Func& func
    )
    {
#ifndef FOAM_USE_TBB
        forAll(list, i)
        {
            func(i);
        }
#else
        try
        {
            tbb::parallel_for
            (
                0,
                list.size(),
                func
            );
        }
        catch (std::exception& ex)
        {
            throw error(ex.what());
        }
#endif
    }
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
