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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::courant

Description
    Reports on maximum and mean Courant number and location of maxima

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    courant.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_courant_H
#define functionObjects_courant_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fields/volFields/volFieldsFwd.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class courant Declaration
\*---------------------------------------------------------------------------*/

class courant
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        word phiName_;

        word rhoName_;

        //- Disallow default bitwise copy construct
        courant(const courant&);

        //- Disallow default bitwise assignment
        void operator=(const courant&);

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

public:

    //- Runtime type information
    TypeName("courant");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        courant
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~courant();


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
