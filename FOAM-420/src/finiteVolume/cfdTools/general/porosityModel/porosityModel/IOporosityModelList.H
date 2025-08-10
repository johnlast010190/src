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
    (c) 2012 OpenFOAM Foundation

Class
    Foam::IOporosityModelList

Description
    List of porosity models with IO functionality

SourceFiles
    IOporosityModelList.C

\*---------------------------------------------------------------------------*/

#ifndef IOporosityModelList_H
#define IOporosityModelList_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IOporosityModelList Declaration
\*---------------------------------------------------------------------------*/

class IOporosityModelList
:
    public IOdictionary,
    public porosityModelList
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh, bool) const;

        //- Disallow default bitwise copy construct
        IOporosityModelList(const IOporosityModelList&);

        //- Disallow default bitwise assignment
        void operator=(const IOporosityModelList&);


public:

    // Constructors

        //- Construct from mesh
        IOporosityModelList
        (
            const fvMesh& mesh,
            bool readFromFvOptions = false
        );

        //- Construct from mesh and registry
        IOporosityModelList
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            bool readFromFvOptions = false
        );


        //- Destructor
        virtual ~IOporosityModelList()
        {}


    // Member Functions

        //- Read dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
