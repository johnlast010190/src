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
    (c) 2015 OpenCFD Ltd.

Class
    Foam::surfaceReader

Description
    Base class for surface readers

SourceFiles
    surfaceReader.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceReader_H
#define surfaceReader_H

#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "MeshedSurface/MeshedSurfaces.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceReader Declaration
\*---------------------------------------------------------------------------*/

class surfaceReader
{
protected:

    //- File name
    fileName fileName_;


public:

    //- Runtime type information
    TypeName("surfaceReader");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceReader,
            fileName,
            (
                const fileName& fName
            ),
            (fName)
        );


    // Selectors

        //- Return a reference to the selected surfaceReader
        static autoPtr<surfaceReader> New
        (
            const word& readType,
            const fileName& fName
        );


    // Constructors

        //- Construct from fileName
        surfaceReader(const fileName& fName);


    //- Destructor
    virtual ~surfaceReader();


    // Member Functions

        //- Return a reference to the surface geometry
        virtual const meshedSurface& geometry() = 0;

        //- Return a list of the available times
        virtual instantList times() const = 0;

        //- Return a list of the available fields at a given time
        virtual wordList fieldNames(const label timeIndex) const = 0;

        //- Return a scalar field at a given time
        virtual tmp<Field<scalar>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const scalar& refValue = pTraits<scalar>::zero
        ) const = 0;

        //- Return a vector field at a given time
        virtual tmp<Field<vector>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const vector& refValue = pTraits<vector>::zero
        ) const = 0;

        //- Return a sphericalTensor field at a given time
        virtual tmp<Field<sphericalTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const sphericalTensor& reValue = pTraits<sphericalTensor>::zero
        ) const = 0;

        //- Return a symmTensor field at a given time
        virtual tmp<Field<symmTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const symmTensor& reValue = pTraits<symmTensor>::zero
        ) const = 0;

        //- Return a tensor field at a given time
        virtual tmp<Field<tensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const tensor& reValue = pTraits<tensor>::zero
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
