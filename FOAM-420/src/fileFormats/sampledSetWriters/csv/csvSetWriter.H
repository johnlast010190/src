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

Class
    Foam::csvSetWriter

Description
    Write set in csv format

SourceFiles
    csvSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef csvSetWriter_H
#define csvSetWriter_H

#include "sampledSetWriters/writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class csvSetWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class csvSetWriter
:
    public writer<Type>
{

    // Private Member Functions

        void writeCoordHeader(const coordSet&, Ostream&) const;

        void writeHeader(const coordSet&, const wordList&, Ostream&) const;


protected:

    virtual void writeSeparator(Ostream&) const override;


public:

    //- Runtime type information
    TypeName("csv");


    // Constructors

        //- Construct null
        csvSetWriter();


    //- Destructor
    virtual ~csvSetWriter();


    // Member Functions

        const char* ext() const override;

        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<const Field<Type>*>&,
            Ostream&
        ) const override;

        virtual void write
        (
            const bool writeTracks,
            const PtrList<coordSet>&,
            const wordList& valueSetNames,
            const List<List<Field<Type>>>&,
            Ostream&
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSetWriters/csv/csvSetWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
