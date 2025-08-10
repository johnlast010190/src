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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/data/data.H"
#include "db/Time/Time.H"
#include "matrices/LduMatrix/LduMatrix/solverPerformance.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::data::setSolverPerformance
(
    const word& name,
    const SolverPerformance<Type>& sp
) const
{
    dictionary& dict = const_cast<dictionary&>(solverPerformanceDict());

    List<SolverPerformance<Type>> perfs;

    bool wordExists = false;
    label wordExistsIter = 0;
    for (const auto& tuple : prevTimeIndices_)
    {
        if (tuple.first() == name)
        {
            wordExists = true;
            break;
        }
        wordExistsIter++;
    }

    if (wordExists)
    {
        if (prevTimeIndices_[wordExistsIter].second()!= this->time().timeIndex())
        {
            prevTimeIndices_[wordExistsIter].second()
                = this->time().timeIndex();
            dict.remove(name);
        }
        else
        {
            dict.readIfPresent(name, perfs);
        }
    }
    else
    {
        prevTimeIndices_.append
            (
                Tuple2<word,label>(name, this->time().timeIndex())
            );
    }

    // Append to list
    perfs.setSize(perfs.size()+1, sp);

    dict.set(name, perfs);
}


template<class Type>
void Foam::data::setSolverPerformance
(
    const SolverPerformance<Type>& sp
) const
{
    setSolverPerformance(sp.fieldName(), sp);
}


template<class Type>
void Foam::data::setBlockSolverPerformance
(
    const word& name,
    const BlockSolverPerformance<Type>& sp
) const
{
    dictionary& dict = const_cast<dictionary&>(blockSolverPerformanceDict());

    List<BlockSolverPerformance<Type>> perfs;

    bool wordExists = false;
    label wordExistsIter = 0;
    for (const auto& tuple : prevTimeBlockIndices_)
    {
        if (tuple.first() == name)
        {
            wordExists = true;
            break;
        }
        wordExistsIter++;
    }

    if (wordExists)
    {
        if
        (
            prevTimeBlockIndices_[wordExistsIter].second()
          !=this->time().timeIndex()
        )
        {
            prevTimeBlockIndices_[wordExistsIter].second()
                = this->time().timeIndex();
            dict.remove(name);
        }
        else
        {
            dict.readIfPresent(name, perfs);
        }
    }
    else
    {
        prevTimeBlockIndices_.append
            (
                Tuple2<word,label>(name, this->time().timeIndex())
            );
    }

    // Append to list
    perfs.setSize(perfs.size()+1, sp);

    dict.set(name, perfs);
}


template<class Type>
void Foam::data::setBlockSolverPerformance
(
    const BlockSolverPerformance<Type>& sp
) const
{
    setBlockSolverPerformance(sp.fieldName(), sp);
}


// ************************************************************************* //
