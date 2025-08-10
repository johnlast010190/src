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
    (c) 2022 Esi Ltd.

Class
    Foam::morphFieldMapper

Description
    Abstract base class to hold the Field mapping for mesh morphs.

\*---------------------------------------------------------------------------*/

#ifndef morphFieldMapper_H
#define morphFieldMapper_H

#include "fields/Fields/fieldMappers/generalFieldMapper/generalFieldMapper.H"
#include "containers/HashTables/Map/Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class morphFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class morphFieldMapper
:
    public generalFieldMapper
{
public:

    // Constructors

        //- Null constructor
        morphFieldMapper()
        {}


    //- Destructor
    virtual ~morphFieldMapper()
    {}


    // Member Functions

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const = 0;

        //- No indirect mapping
        virtual bool indirect() const override
        {
            return false;
        }

        //- Are there any inserted objects
        virtual bool insertedObjects() const = 0;

        //- Return list of inserted objects
        virtual const labelList& insertedObjectLabels() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
