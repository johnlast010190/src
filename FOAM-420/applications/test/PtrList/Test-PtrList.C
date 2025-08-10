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

Application

Description

\*---------------------------------------------------------------------------*/

#include "include/OSspecific.H"

#include "primitives/Scalar/scalar/scalar.H"
#include "db/IOstreams/IOstreams.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "meshes/primitiveShapes/plane/plane.H"
#include "containers/Lists/DynamicList/DynamicList.H"

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

    autoPtr<Scalar> clone() const
    {
        return autoPtr<Scalar>(new Scalar(data_));
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
    PtrList<Scalar> list1(10);
    PtrList<Scalar> list2(15);
    PtrList<Scalar> listApp;

    forAll(list1, i)
    {
        list1.set(i, new Scalar(1.3*i));
    }

    forAll(list2, i)
    {
        list2.set(i, new Scalar(10 + 1.3*i));
    }

    for (label i = 0; i < 5; ++i)
    {
        listApp.append(new Scalar(1.3*i));
    }

    Info<<"list1: " << list1 << endl;
    Info<<"list2: " << list2 << endl;
    Info<<"listApp: " << listApp << endl;

    Info<<"indirectly delete some items via set(.., 0) :" << endl;
    for (label i = 0; i < 3; i++)
    {
        list1.set(i, 0);
    }

    Info<<"transfer list2 -> list1:" << endl;
    list1.transfer(list2);

    Info<<"list1: " << list1 << endl;
    Info<<"list2: " << list2 << endl;

    Info<<"indirectly delete some items via setSize :" << endl;
    list1.setSize(4);

    Info<<"list1: " << list1 << endl;

    PtrList<Scalar> list3(list1.xfer());
    Info<< "Transferred via the xfer() method" << endl;

    Info<<"list1: " << list1 << endl;
    Info<<"list2: " << list2 << endl;
    Info<<"list3: " << list3 << endl;

    PtrList<plane> planes;
    planes.append(new plane(vector::one, vector::one));
    planes.append(new plane(vector(1,2,3), vector::one));

    forAll(planes, p)
        Info<< "plane " << planes[p] << endl;

    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
