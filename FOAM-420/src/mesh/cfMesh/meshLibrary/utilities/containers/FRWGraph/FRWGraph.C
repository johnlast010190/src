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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

#include "utilities/containers/FRWGraph/FRWGraph.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, Foam::label width>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::FRWGraph<T, width>& DL
)
{
    os << DL.size() << "(" << nl;

    for (label i=0;i<DL.size();++i)
    {
        os << width << "(";

        for (label j=0;j<width;++j)
        {
            if (j)
            {
                os << " ";
            }

            os << DL(i, j);
        }

        os << ")" << nl;
    }

    os << ")";

    // Check state of IOstream
    os.check
    (
        "template<class T, Foam::label width>Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream& os, const Foam::FRWGraph<T, width>&)"
    );

    return os;
}

/*
template<class T, Foam::label width>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::FRWGraph<T, width>& DL
)
{
    label size;
    T e;
    is >> size;
    DL.setSize(size);
    for (IndexType i=0;i<size;++i)
    {
        is >> e;
        DL[i] = e;
    }

    return is;
}
*/

// ************************************************************************* //
