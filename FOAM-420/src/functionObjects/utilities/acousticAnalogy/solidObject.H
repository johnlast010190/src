
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
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::kfwhModel

Group
    grpFunctionObjects

Description
	Define primitive solid object for identifying whether part of a surface
	patch is inside a solid object

SourceFiles
    solidObject.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/IOstreams/Fstreams/OFstream.H"
#include "primitives/Vector/vector/vector.H"
#include "fvMesh/fvMesh.H"

#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "smBoundBox.H"

namespace Foam
{

class solidObject
{
	public:
		word name_;
		word type_;
		scalar r_;
		vectorField points_;
        autoPtr<smBoundBox> solidBox_;
	public:
		//- Runtime type information
		TypeName("solidObject");

		declareRunTimeSelectionTable
        (
            autoPtr,
            solidObject,
            word,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );

        static autoPtr<solidObject> New
        (
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);

		inline bool inObject(const point &pt);

		solidObject();

		solidObject
		(
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);

        const word& name() const
        {
            return name_;
        }

		virtual ~solidObject(){;}
};

bool solidObject::inObject(const point &p)
{
	if (type_=="cylinder")
	{
		if (points_.size()<2)
		{
			return false;
		}

		const point &a=points_[0];
		const point &b=points_[1];
		point n=(b-a)/(mag(b-a)+1.0e-20);
		scalar d1=(p-b)&n;
		if (d1>0)
		{
			return false;
		}
		scalar d2=(p-a)&n;
		if (d2<0)
		{
			return false;
		}
	    point pt=(a-p)-((a-p)&n)*n;
	    scalar d=mag(pt);
	    if (d<r_)
	    {
			return true;
		}
		return false;
	}

	if (type_=="sphere")
	{
		const point &a=points_[0];
		scalar r1=mag(p-a);
		if (r1<r_)
		{
			return true;
		}
		return false;
	}

    if (type_=="smBoundBox")
	{
        if (solidBox_.empty())
        {
            return false;
        }
        if (!solidBox_().contains(p))
        {
            return false;
        }
        if (solidBox_().pointInBox(p))
        {
            return false; //in a fluid grid
        }
        else
        {
            return true; //in a solid grid
        }
	}

	if (type_=="box")
	{
		if (points_.size()<2)
		{
			FatalErrorInFunction
            << "Number of points must be at least 2 to define a box"
            << exit(FatalError);
		}
		const point &a=points_[0];
		const point &b=points_[1];
		boundBox bbox(a,b);
		if (bbox.contains(p))
		{
			return true;
		}
		return false;
	}

	return false;
}


}//Foam
