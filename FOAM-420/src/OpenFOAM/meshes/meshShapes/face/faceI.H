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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::face::face()
{}


inline Foam::face::face(const label sz)
:
    labelList(sz, -1)
{}


inline Foam::face::face(const labelUList& lst)
:
    labelList(lst)
{}


inline Foam::face::face(std::initializer_list<label> lst)
:
    labelList(lst)
{}


inline Foam::face::face(const Xfer<labelList>& lst)
:
    labelList(lst)
{}


inline Foam::face::face(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::pointField Foam::face::points
(
    const UList<point>& meshPoints
) const
{
    // There are as many points as there are labels for them
    pointField p(size());

    // For each point in list, set it to the point in 'pnts' addressed
    // by 'labs'
    forAll(p, i)
    {
        p[i] = meshPoints[operator[](i)];
    }

    // Return list
    return p;
}


inline Foam::vector Foam::face::unitNormal(const UList<point>& p) const
{
    const vector n(areaNormal(p));
    const scalar s(Foam::mag(n));
    return s < ROOTVSMALL ? Zero : n/s;
}


inline Foam::scalar Foam::face::mag(const UList<point>& p) const
{
    return ::Foam::mag(areaNormal(p));
}


inline Foam::label Foam::face::nEdges() const
{
    // for a closed polygon a number of edges is the same as number of points
    return size();
}


inline Foam::edge Foam::face::faceEdge(const label n) const
{
    return edge(operator[](n), operator[](fcIndex(n)));
}


inline bool Foam::face::found(const label globalIndex) const
{
    return which(globalIndex) != -1;
}


inline Foam::label Foam::face::nextLabel(const label i) const
{
    return operator[](fcIndex(i));
}


inline Foam::label Foam::face::prevLabel(const label i) const
{
    return operator[](rcIndex(i));
}


inline Foam::label Foam::face::nTriangles() const
{
    return size() - 2;
}



// * * * * * * * * * * * * * * * Global Operators   * * * * * * * * * * * * * //

inline bool Foam::operator==(const face& a, const face& b)
{
    return face::compare(a,b) != 0;
}

inline bool Foam::operator!=(const face& a, const face& b)
{
    return face::compare(a,b) == 0;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, face& f)
{
    is >> static_cast<labelList&>(f);
    is.check(FUNCTION_NAME);
    return is;
}

// ************************************************************************* //
