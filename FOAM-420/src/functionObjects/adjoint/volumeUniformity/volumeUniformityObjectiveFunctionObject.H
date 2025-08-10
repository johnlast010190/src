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
    volumeUniformityObjectiveFunctionObject

Description
    Calculates the flow uniformity inside a volume.

SourceFiles
    volumeUniformityObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef volumeUniformityObjectiveFunctionObject_H
#define volumeUniformityObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class volumeUniformityObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class volumeUniformityObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- Uniformity in the secified direction
        autoPtr<vector> desiredFlowDirection_;

        //- Specified target velocity
        autoPtr<vector> targetVelocity_;

        //- Objective multiplier
        scalar sourceCoeff_;

        //- Zone volume
        scalar volume_;


    // Protected Member Functions

        scalar zoneVolume() const;

        vector targetVelocity() const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        volumeUniformityObjectiveFunctionObject
        (
            const volumeUniformityObjectiveFunctionObject&
        );
        void operator=(const volumeUniformityObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("volumeUniformity");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        volumeUniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        volumeUniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~volumeUniformityObjectiveFunctionObject();


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
