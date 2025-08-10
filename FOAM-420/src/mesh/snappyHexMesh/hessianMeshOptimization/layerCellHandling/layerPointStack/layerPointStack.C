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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/layerCellHandling/layerPointStack/layerPointStack.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    layerPointStack::layerPointStack
    (
        const bool& baseStack,
        const List<vector>& stackPoints,
        const labelList& pointStackLabels,
        const label& layerNu
    )
    :
        baseStack_(baseStack),
        stackPoints_(stackPoints),
        pointStackLabels_(pointStackLabels),
        syncAddressing_(-1),
        layerNu_(layerNu)
    {
    };

    layerPointStack::layerPointStack()
    :
        baseStack_(false),
        stackPoints_(0),
        pointStackLabels_(0),
        syncAddressing_(-1),
        layerNu_(0)
    {
    };

    void layerPointStack::setAddress(const label& add)
    {
        syncAddressing_ = add;
    }

    void layerPointStack::increment(const vector& p)
    {
        stackPoints_.append(p);
    }

    void layerPointStack::substitute(const List<vector>& newStack)
    {
        stackPoints_ = newStack;
    }

    void layerPointStack::merge(const List<layerPointStack>& lPSList)
    {
        forAll(pointStackLabels_, pI)
        {
            const label& pL = pointStackLabels_[pI];
            const layerPointStack& lPS = lPSList[pL];
            if (this!=&lPS)
            {
                merge(*this, lPS);
            }
        }
    }

    void layerPointStack::merge(layerPointStack& ps1, const layerPointStack& ps2)
    {
        const List<vector>& stack1 = ps1.stackPoints();
        const List<vector>& stack2 = ps2.stackPoints();
        label start = -1;
        label end = -1;
        if (stack1.size())
        {
            forAll(stack2, pI)
            {
                if (stack1[0] == stack2[pI])
                {
                    start = pI;
                }

                if (stack1[stack1.size()-1] == stack2[pI])
                {
                    end = pI;
                }
            }
        }
        else
        {
            ps1.substitute(ps2.stackPoints());
        }
        if (start != -1 && end != -1)
        {
            ps1.substitute(ps2.stackPoints());
        }
        else if (start > 0 && end == -1)
        {
            List<vector> newStackPoint(ps1.size()+start);
            for (int i = 0; i<start; i++)
            {
                newStackPoint[i] = ps2.stackPoints()[i];
            }
            for (int i = start; i<newStackPoint.size(); i++)
            {
                int j = i - start;
                newStackPoint[i] = ps1.stackPoints()[j];
            }
            ps1.substitute(newStackPoint);
        }
        else if (start == -1 && end > -1)
        {
            List<vector> newStackPoint(ps1.size()+ps2.size()-end-1);
            for (int i = 0; i<ps1.size(); i++)
            {
                newStackPoint[i] = ps1.stackPoints()[i];
            }
            for (int i = ps1.size(); i<newStackPoint.size(); i++)
            {
                int j = i - ps1.size() + end + 1;
                newStackPoint[i] = ps2.stackPoints()[j];
            }
            ps1.substitute(newStackPoint);
        }
    }

    label layerPointStack::getIndex(const vector& point) const
    {
        forAll(stackPoints(), pI)
        {
            if
            (
                mag(point.x()-stackPoints()[pI].x())<1.e-10
             &&    mag(point.y()-stackPoints()[pI].y())<1.e-10
             && mag(point.z()-stackPoints()[pI].z())<1.e-10
            )
            {
                return pI;
            }
        }
        return -1;
    }
};
