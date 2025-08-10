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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::faceSelections::faceZoneSelection

Description
    Select faces from faceZone

SourceFiles
    faceZoneSelection.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneSelection_H
#define faceZoneSelection_H

#include "faceSelection.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/bools/lists/boolList.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionSplit;

namespace faceSelections
{

/*---------------------------------------------------------------------------*\
                     Class faceZoneSelection Declaration
\*---------------------------------------------------------------------------*/

class faceZoneSelection
:
    public faceSelection
{
    // Private data

        //- Name of faceZone
        const word zoneName_;


   // Private Member Functions

public:

    //- Runtime type information
    TypeName("faceZone");


    // Constructors

        //- Construct from dictionary
        faceZoneSelection
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Clone
        autoPtr<faceSelection> clone() const
        {
            NotImplemented;
        }


    //- Destructor
    virtual ~faceZoneSelection();


    // Member Functions

        //- Apply this selector
        virtual void select(const label zoneID, labelList&, boolList&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace faceSelections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
