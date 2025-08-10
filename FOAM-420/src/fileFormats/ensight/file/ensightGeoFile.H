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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::ensightGeoFile

Description
    Specialized Ensight output with extra geometry file header

\*---------------------------------------------------------------------------*/

#ifndef ensightGeoFile_H
#define ensightGeoFile_H

#include "ensight/file/ensightFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightGeoFile Declaration
\*---------------------------------------------------------------------------*/

class ensightGeoFile
:
    public ensightFile
{
    // Private Member Functions

        //- Initialize by outputting header information
        void initialize();

        //- Disallow default bitwise assignment
        void operator=(const ensightGeoFile&) = delete;

        //- Disallow default copy constructor
        ensightGeoFile(const ensightGeoFile&) = delete;


public:

    // Constructors

        //- Construct from pathname.
        //  The entire pathname is checked for valid ensight naming.
        ensightGeoFile
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::BINARY
        );

        //- Construct from path and name.
        //  Only the name portion is checked for valid ensight naming.
        ensightGeoFile
        (
            const fileName& path,
            const fileName& name,
            IOstream::streamFormat format=IOstream::BINARY
        );


    //- Destructor
    ~ensightGeoFile();


    // Output

        //- Write keyword with trailing newline
        virtual Ostream& writeKeyword(const keyType& key);


    // Convenience Output Methods

        using ensightFile::beginPart;

        //- Begin a "part" (0-based index), with a description.
        void beginPart(const label index, const string& description);

        //- Begin a "coordinates" block
        void beginCoordinates(const label npoints);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
