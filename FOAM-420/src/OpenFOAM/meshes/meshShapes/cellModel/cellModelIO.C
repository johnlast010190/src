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

#include "meshes/meshShapes/cellModel/cellModel.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::cellModel::cellModel(Istream& is)
{
    dictionaryEntry entry(dictionary::null, is);
    name_ = entry.keyword();
    entry.lookup("index") >> index_;
    entry.lookup("numberOfPoints") >> nPoints_;
    entry.lookup("faces") >> faces_;
    entry.lookup("edges") >> edges_;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const cellModel& c)
{
    os  << "name" << tab << c.name_ << tab
        << "index" << tab << c.index_ << tab
        << "numberOfPoints" << tab << c.nPoints_ << tab
        << "faces" << tab << c.faces_ << tab
        << "edges" << tab << c.edges_ << endl;

    return os;
}


template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<cellModel>& ip)
{
    const cellModel& cm = ip.t_;

    os  << "name = " << cm.name() << ", "
        << "index = " << cm.index() << ", "
        << "number of points = " << cm.nPoints() << ", "
        << "number of faces = " << cm.nFaces() << ", "
        << "number of edges = " << cm.nEdges()
        << endl;

    return os;
}


// ************************************************************************* //
