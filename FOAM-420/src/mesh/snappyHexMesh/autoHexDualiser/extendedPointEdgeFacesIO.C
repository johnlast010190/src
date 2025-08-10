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

#include "autoHexDualiser/extendedPointEdgeFaces.C"
#include "db/IOstreams/IOstreams/Istream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/token/token.H"
#include "primitives/contiguous/contiguous.H"


namespace Foam
{
    // Stream operators

    extendedPointEdgeFaces::extendedPointEdgeFaces(Istream& is)
    {
        operator>>(is, *this);
    }

    Istream& operator>>(Istream& is, extendedPointEdgeFaces& pI)
    {
        faceList faces;
        pointList pts;
        labelList addr;
        wordList names;
        labelList neighbours;

        is >> faces >> pts >> addr >> names >> neighbours;

        pI = extendedPointEdgeFaces(faces,pts,addr,names,neighbours);

        return is;
    }

    Ostream& operator<<(Ostream& os, const extendedPointEdgeFaces& pI)
    {
        os << pI.edgeFaces() <<" "<< pI.points()<<" "<<pI.addressing()
           <<" "<<pI.edgePatchNames()<<" "<<pI.neighbours();

        return os;
    }
};
