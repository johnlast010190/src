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
    (c) 2009-2021 Esi Ltd.

Class
    Foam::dataItem

Description
    A database item object for machine-learning

SourceFiles
    dataItem.C

\*---------------------------------------------------------------------------*/


#pragma once

//#include <unistd.h>
#include <algorithm>
#include <sstream>
#include "foamMap.H"
#include "array3d.H"
#include "gridField.H"
#include "blockSample.H"

namespace Foam
{

inline void split
(
	std::vector<word>& vect,
	const word& str
)
{
	vect.resize(0);
    std::istringstream iss(str);
    for (word w; iss >> w;)
    {
        vect.push_back(w);
    }
}

inline void split
(
	std::vector<word>& elems,
	const word& s,
	char delim
)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
}

class dataItem
{

protected:
    word caseName_;
    scalar rhoRef_;
    vector Uref_;
    word caseType_;
    word caseDir_;

    scalar identity_; //label of the case as training class label
    scalar xmin_;
    scalar ymin_;
    scalar zmin_;
    scalar xmax_;
    scalar ymax_;
    scalar zmax_;
    label m_;
    label n_;
    label k_;

public:

    void setDataInfo(const word &dat);
    std::vector<scalar> geoData_;

    label m() const
    {
        return m_;
    }

    label n() const
    {
        return n_;
    }

    label k() const
    {
        return k_;
    }

    dataItem(const word& dat);

    void setData
    (
        label m,
        label n,
        label k,
        const std::vector<scalar>& data
    );

    scalar getDist
    (
        const std::vector<scalar>& vdata
    ) const;

    scalar gdataSum() const;

    word path() const
    {
        return caseDir_;
    }

    void setData
    (
        const word& data
    );

    label numSolids() const;

    scalar identity() const
    {
        return identity_;
    }

    const word &caseName() const
    {
        return caseName_;
    }

    void identity(scalar d)
    {
        identity_=d;
    }

    void show() const;

    void clearGeoData()
    {
        geoData_.resize(0);
    }

    dataItem ();
    ~dataItem (){;}
};


}//Foam
