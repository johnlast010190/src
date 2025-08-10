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
    (c) 2017-2021 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::cylinderToFace

Description
    A topoSetSource to select faces based on face centres inside a cylinder.

SourceFiles
    cylinderToFace.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderToFace_H
#define cylinderToFace_H

#include "sets/topoSetSource/topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cylinderToFace Declaration
\*---------------------------------------------------------------------------*/

class cylinderToFace
:
    public topoSetSource
{
    // Private Data

        //- First point on cylinder axis
        vector point1_;

        //- Second point on cylinder axis
        vector point2_;

        //- Radius
        scalar radius_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("cylinderToFace");


    // Constructors

        //- Construct from components
        cylinderToFace
        (
            const polyMesh& mesh,
            const vector& point1,
            const vector& point2,
            const scalar radius
        );

        //- Construct from dictionary
        cylinderToFace(const polyMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~cylinderToFace();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACESETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
