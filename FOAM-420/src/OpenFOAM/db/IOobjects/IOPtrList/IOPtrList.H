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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2018 OpenCFD Ltd.

Class
    Foam::IOPtrList

Description
    A PtrList of objects of type \<T\> with automated input and output.

SourceFiles
    IOPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef IOPtrList_H
#define IOPtrList_H

#include "containers/Lists/PtrList/PtrList.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOPtrList
:
    public regIOobject,
    public PtrList<T>
{

public:

    //- Runtime type information
    TypeName("PtrList");


    // Constructors

        //- Default copy construct
        IOPtrList(const IOPtrList&) = default;

        //- Construct from IOobject using given Istream constructor class
        template<class INew>
        IOPtrList(const IOobject&, const INew&);

        //- Construct from IOobject
        explicit IOPtrList(const IOobject&);

        //- Construct from IOobject with given size
        IOPtrList(const IOobject&, const label);

        //- Construct from IOobject and a PtrList
        IOPtrList(const IOobject&, const PtrList<T>&);

        //- Construct by transferring the PtrList contents
        IOPtrList(const IOobject&, const Xfer<PtrList<T>>&);


    //- Destructor
    virtual ~IOPtrList();


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        //- Copy assignment of entries
        void operator=(const IOPtrList<T>&);

        //- Copy or move assignment of entries
        using PtrList<T>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/IOobjects/IOPtrList/IOPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
