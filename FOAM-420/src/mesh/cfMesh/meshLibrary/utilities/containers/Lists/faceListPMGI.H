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

Description


\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/Lists/SubList/SubList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Constructors

inline faceListPMG::faceListPMG(const IOobject& io)
:
    regIOobject(io),
    faceList(readStream(typeName)),
    nElmts_(faceList::size())
{
}

inline faceListPMG::faceListPMG(const IOobject& io, const label s)
:
    regIOobject(io),
    faceList(s),
    nElmts_(s)
{
}

inline faceListPMG::faceListPMG(const IOobject& io, const faceList& fcs)
:
    regIOobject(io),
    faceList(fcs),
    nElmts_(fcs.size())
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor

inline faceListPMG::~faceListPMG()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label faceListPMG::size() const
{
    return nElmts_;
}

inline void faceListPMG::setSize(const label nElmts)
{
    if (nElmts >= faceList::size())
    {
        if (faceList::size() != 0)
        {
            Info<< "Resizing faces!" << endl;
            faceList copy(label(1.5*nElmts));
            for (label i=0;i<nElmts_;++i)
                copy[i].transfer(this->operator[](i));

            faceList::transfer(copy);
        }
        else
        {
            faceList::setSize(label(1.5*nElmts));
        }
    }

    nElmts_ = nElmts;
}

inline void faceListPMG::clear()
{
    nElmts_ = 0;
}

inline void faceListPMG::append(const face& f)
{
    const label i = nElmts_;
    setSize(i+1);
    this->operator[](i) = f;
}

inline face& faceListPMG::newElmt(const label fI)
{
    setSize(fI+1);
    return this->operator[](fI);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool faceListPMG::writeData(Ostream& os) const
{
    return (os << *this).good();
}

inline void faceListPMG::operator=(const faceList& fcs)
{
    setSize(fcs.size());
    forAll(fcs, fI)
        this->operator[](fI) = fcs[fI];
}

inline Ostream& operator<<(Ostream& os, const faceListPMG& fcs)
{
    SubList<face> f(fcs, fcs.nElmts_, 0);

    os << f;
    return os;
}

inline Istream& operator>>(Istream& is, faceListPMG& fcs)
{
    faceList& faces = static_cast<faceList&>(fcs);
    is >> faces;
    fcs.nElmts_  = faces.size();

    return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
