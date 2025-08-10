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
    Foam::setUpdater

Description
    Keeps cell/face/vertex sets uptodate. Both the ones loaded and the ones
    on disk.

SourceFiles
    setUpdater.C
    setUpdaterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef setUpdater_H
#define setUpdater_H

#include "polyTopoChange/polyMeshModifier/polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class setUpdater Declaration
\*---------------------------------------------------------------------------*/

class setUpdater
:
    public polyMeshModifier
{
    // Private Member Functions

        //- Updates all sets
        template<class Type>
        void updateSets(const mapPolyMesh& morphMap) const;

        //- Disallow default bitwise copy construct
        setUpdater(const setUpdater&);

        //- Disallow default bitwise assignment
        void operator=(const setUpdater&);


public:

    //- Runtime type information
    TypeName("setUpdater");


    // Constructors

        //- Construct from dictionary
        setUpdater
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    //- Destructor
    virtual ~setUpdater();


    // Member Functions

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "setUpdater/setUpdaterTemplates.C"
#endif


#endif

// ************************************************************************* //
