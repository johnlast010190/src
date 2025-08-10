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
    Foam::f4stSetWriter

SourceFiles
    f4stSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef f4stSetWriter_H
#define f4stSetWriter_H

#include "sampledSetWriters/writer.H"
#include "primitives/Vector/vector/vector.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/IOstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class f4stSetWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class f4stSetWriter: public writer<Type>
{

// Private Member Functions

protected:

    bool binary_;

    bool compressed_;

    mutable bool meshIsWritten_;

    mutable autoPtr<OFstream> outFilePtr_;

    mutable DynamicList<scalar> timeList_;

    void createFile() const;

    void closeFile() const;

    void writePoints(const coordSet& points) const;

    void writeData
    (
        const wordList& valueSetNames,
        const List<const Field<Type>*>& valueSets,
        const scalar& time
    ) const;

    void writeTimeValue(const scalar& timeValue) const;

    bool meshIsWritten() const;

    bool timeIsWritten(const scalar& timeValue) const;

public:

//- Runtime type information
TypeName("f4st");


// Constructors

//- Construct null
f4stSetWriter();


//- Destructor
virtual ~f4stSetWriter();


// Member Functions

    const char* ext() const override;

    virtual bool sepFile() const override {return true;};

    virtual void write
    (
        const coordSet&,
        const wordList&,
        const List<const Field<Type>*>&,
        const scalar
    ) const override;

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
    ) const override
{
    NotImplemented;
}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "sampledSetWriters/f4st/f4stSetWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
