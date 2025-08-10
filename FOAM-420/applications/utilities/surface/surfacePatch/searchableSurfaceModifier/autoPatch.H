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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::searchableSurfaceModifier::autoPatch

Description
    Patchify triangles based on a feature angle. Replacement of surfaceAutoPatch

SourceFiles
    autoPatch.C

\*---------------------------------------------------------------------------*/

#ifndef autoPatch_H
#define autoPatch_H

#include "searchableSurfaceModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace searchableSurfaceModifiers
{

/*---------------------------------------------------------------------------*\
                          Class autoPatch Declaration
\*---------------------------------------------------------------------------*/

class autoPatch
:
    public searchableSurfaceModifier
{
    // Private data

        //- Feature angle
        const scalar featureAngle_;


   // Private Member Functions

public:

    //- Runtime type information
    TypeName("autoPatch");


    // Constructors

        //- Construct from dictionary
        autoPatch(const searchableSurfaces&, const dictionary&);

        //- Clone
        autoPtr<searchableSurfaceModifier> clone() const
        {
            NotImplemented;
        }


    //- Destructor
    virtual ~autoPatch();


    // Member Functions

        //- Apply this selector
        virtual bool modify(const labelList& regions, searchableSurface&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace searchableSurfaceModifiers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
