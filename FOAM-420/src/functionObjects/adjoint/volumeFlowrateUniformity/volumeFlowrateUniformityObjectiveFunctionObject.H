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
    (c) 2019 Esi Ltd.

Authors
    2019. Nikolaos Magoulas (Esi Ltd.). All rights reserved.

Class
    Foam::functionObjects::objectiveFunctionObject::
    volumeFlowrateUniformityObjectiveFunctionObject

Description
    Calculates the volume flow rate split in the specified direction.

SourceFiles
    volumeFlowrateUniformityObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef volumeFlowrateUniformityObjectiveFunctionObject_H
#define volumeFlowrateUniformityObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class volumeFlowrateUniformityObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class volumeFlowrateUniformityObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- Mass flow uniformity in the specified direction
        vectorField flowDirection_;

        //- Zone maximum thickness in the flow direction
        scalarList ts_;

        //- Zone frontal flow area
        scalarList areas_;

        //- Zone target mass flow fraction
        scalarList zoneTargetFraction_;


private:

    // Private Member Functions

        //- Calculate single zone's volume
        scalar zoneVolume(word) const;

        //- Calculate zone's thickness in the vector direction
        scalar zoneThickness(word, vector) const;

        //- Disallow default bitwise copy construct and assignment
        volumeFlowrateUniformityObjectiveFunctionObject
        (
            const volumeFlowrateUniformityObjectiveFunctionObject&
        );
        void operator=(const volumeFlowrateUniformityObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("volumeFlowrateUniformity");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        volumeFlowrateUniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        volumeFlowrateUniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~volumeFlowrateUniformityObjectiveFunctionObject();


    // Member Functions

        //- Read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate the value of the cost function and write it to file
        virtual bool execute();

        //- Track and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
