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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldInstance/fieldInstance.H"

namespace Foam
{
    // Stream operators

    Istream& operator>> (Istream& is, fieldInstance& fI)
    {
        FOAM_ASSERT(is.format() != IOstream::ASCII)
        {
            FatalErrorInFunction << "ASCII fieldInstance input is not supported." << abort(FatalError);
        }

        scalar v(readScalar(is));
        point p(is);
        fI = fieldInstance(v,p);
        return is;
    }

    Ostream& operator<< (Ostream& os, const fieldInstance& fI)
    {
        // Text-mode seems to be used for writing this stuff to the terminal,
        // so regrettably we have to stuff this annoying branch in here.
        if (os.format() == IOstream::ASCII) {
            os << fI.value() << " " << fI.position();
        } else {
            os << fI.value() << fI.position();
        }

        return os;
    }

};


// ************************************************************************* //
