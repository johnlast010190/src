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
    (c) 2016-2018 OpenCFD Ltd.
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::GlobalIOField

Description
    IOField with global data (so optionally read from master)

SourceFiles
    GlobalIOField.C

\*---------------------------------------------------------------------------*/

#ifndef GlobalIOField_H
#define GlobalIOField_H

#include "db/regIOobject/regIOobject.H"
#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GlobalIOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GlobalIOField
:
    public regIOobject,
    public Field<Type>
{

public:

    TypeName("Field");


    // Constructors

        //- Default copy construct
        GlobalIOField(const GlobalIOField&) = default;

        //- Construct from IOobject
        GlobalIOField(const IOobject&);

        //- Construct from IOobject and size (does not set values)
        GlobalIOField(const IOobject&, const label size);

        //- Construct from components
        GlobalIOField(const IOobject&, const Field<Type>&);

        //- Construct by transferring the Field contents
        GlobalIOField(const IOobject&, const Xfer<Field<Type>>&);


    //- Destructor
    virtual ~GlobalIOField();


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const GlobalIOField<Type>&);

        void operator=(const Field<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "db/IOobjects/GlobalIOField/GlobalIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
