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
    Foam::writer

Description
    Base class for graphics format writing. Entry points are
    - write(..). \n
      Write to an Ostream a table of points with corresponding values.
    - write(scalar/vector/sphericalTensor/symmTensor/tensor). \n
      Write single scalar/vector/sphericalTensor/symmTensor/tensor.
      Default is to write space separated components.

    Example:
    \verbatim
        // Construct writer of xmgr type
        autoPtr<writer<scalar>> scalarFormatter(writer<scalar>::New("xmgr"));

        // Output list of points and corresponding values
        scalarFormatter().write
        (
            coordSet(...)
            "U.component(0)",   // name of values
            vals                // values
        );
    \endverbatim

SourceFiles
    writer.C

\*---------------------------------------------------------------------------*/

#ifndef writer_H
#define writer_H

#include "primitives/strings/fileName/fileName.H"
#include "primitives/strings/lists/wordList.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/autoPtr/autoPtr.H"
#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class coordSet;

/*---------------------------------------------------------------------------*\
                           Class writer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class writer
{

protected:

    //- Generates filename from coordSet and sampled fields
    fileName getBaseName(const word& ptName, const wordList&) const;

    void writeCoord(const coordSet&, const label sampleI, Ostream&) const;

    //- Writes single-column ascii write. Column 1 is coordSet coordinate,
    //  columns 2 is the value. Uses write() function
    //  to write coordinate in correct format.
    void writeTable(const coordSet&, const List<Type>&, Ostream&) const;

    //- Writes multi-column ascii write. Column 1 is coordSet coordinate,
    //  columns 2..n are the values. Uses write() function
    //  to write coordinate in correct format.
    void writeTable
    (
        const coordSet&,
        const List<const List<Type>*>&,
        Ostream& os
    ) const;

    //- Writes a separator. Used by write functions.
    virtual void writeSeparator(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("writer");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            writer,
            word,
            (),
            ()
        );


    // Selectors

        //- Return a reference to the selected writer
        static autoPtr<writer> New(const word& writeFormat);


    // Constructors

        //- Construct null
        writer();


    //- Destructor
    virtual ~writer() = 0;


    // Member Functions

        //- Generate file name with correct extension
        virtual bool sepFile()  const {return false;}

        virtual const char* ext() const = 0;

        //- Generate file name with correct extension
        fileName getFileName
        (
            const coordSet& points,
            const wordList& valueSetNames
        ) const;

        fileName getFileName
        (
            const word& baseName,
            const wordList& valueSetNames
        ) const;

        //- General entry point for writing.
        //  The data is organized in a set of point with one or more values
        //  per point
        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<const Field<Type>*>&,
            const scalar
        ) const {
                // Do nothing
                } ;

        //- General entry point for writing.
        //  The data is organized in a set of point with one or more values
        //  per point
        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<const Field<Type>*>&,
            Ostream&
        ) const = 0;

        //- General entry point for writing.
        //  The data is organized in a set of point with one or more values
        //  per point
        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<Field<Type>>&,
            Ostream&
        ) const;

        //- General entry point for writing of multiple coordSets.
        //  Each coordSet (track) has same data variables.
        //  The data is per variable, per track, per point of track.
        //  If writeTracks adds connecting lines (wherever applicable)
        virtual void write
        (
            const bool writeTracks,
            const PtrList<coordSet>&,
            const wordList& valueSetNames,
            const List<List<Field<Type>>>&,
            Ostream&
        ) const = 0;

        //- Write scalar as ascii
        virtual Ostream& write(const scalar, Ostream&) const;

        template<class VSType>
        Ostream& writeVS(const VSType&, Ostream&) const;

        //- Write vector. Tab separated ascii
        virtual Ostream& write(const vector&, Ostream&) const;

        //- Write sphericalTensor. Tab separated ascii
        virtual Ostream& write(const sphericalTensor&, Ostream&) const;

        //- Write symmTensor. Tab separated ascii
        virtual Ostream& write(const symmTensor&, Ostream&) const;

        //- Write tensor. Tab separated ascii
        virtual Ostream& write(const tensor&, Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSetWriters/writer.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
