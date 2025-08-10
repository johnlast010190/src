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
    (c) 2015 OpenFOAM Foundation

Application

Description

\*---------------------------------------------------------------------------*/

#include "include/OSspecific.H"

#include "primitives/Scalar/scalar/scalar.H"

#include "db/IOstreams/IOstreams.H"
#include "containers/Dictionaries/PtrListDictionary/PtrListDictionary.H"

using namespace Foam;

class Scalar
{
    scalar data_;

public:

    Scalar()
    :
        data_(0)
    {}

    Scalar(scalar val)
    :
        data_(val)
    {}

    ~Scalar()
    {
        Info<<"delete Scalar: " << data_ << endl;
    }

    friend Ostream& operator<<(Ostream& os, const Scalar& val)
    {
        os  << val.data_;
        return os;
    }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    PtrListDictionary<Scalar> scalarDict(10);
    forAll(scalarDict, i)
    {
        word key("ent" + name(i));
        scalarDict.set(i, key, new Scalar(1.3*i));
    }

    Info<< nl << "scalarDict1: " << endl;
    forAll(scalarDict, i)
    {
        Info<< "elem " << i << " = " << scalarDict[i] << endl;
    }

    Scalar* ent8Ptr = scalarDict.lookupPtr("ent8");

    Info<< "ent8 = " << *ent8Ptr << endl;

    PtrListDictionary<Scalar> scalarDict2(15);
    forAll(scalarDict2, i)
    {
        word key("ent" + name(i));
        scalarDict2.set(i, key, new Scalar(1.3*i));
    }
    Info<< nl << "scalarDict2: " << endl;
    forAll(scalarDict2, i)
    {
        Info<< "elem " << i << " = " << scalarDict2[i] << endl;
    }

    scalarDict.transfer(scalarDict2);

    Scalar* p = scalarDict.lookupPtr("ent8");

    if (p)
    {
        Info<< "found: " << *p << endl;
    }
    else
    {
        Info<< "no p: " << endl;
    }

    scalarDict.clear();

    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
