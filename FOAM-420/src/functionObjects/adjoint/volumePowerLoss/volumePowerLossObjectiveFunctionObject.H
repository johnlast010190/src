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
    volumePowerLossObjectiveFunctionObject

Description
    Calculate pressure loss defined on the volume.

SourceFiles
    volumePowerLossObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef volumePowerLossObjectiveFunctionObject_H
#define volumePowerLossObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class volumePowerLossObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class volumePowerLossObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- total objective patch inlet area
        scalar AreaObjPatchInlet_;

        //- total objective patch outlet area
        scalar AreaObjPatchOutlet_;

        //- Switch on adiabatic reporting for compressible flow (default false)
        Switch adiabatic_;

        //- Switch on tangential boundary stress energy transport
        // (default false)
        Switch stressWork_;

        Switch printLossesSeparately_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        volumePowerLossObjectiveFunctionObject
        (
            const volumePowerLossObjectiveFunctionObject&
        );
        void operator=(const volumePowerLossObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("volumePowerLoss");


    // Constructors

        //- Construct from time, dictionary, adjoint solver and use old-style
        //  adjoint file-format
        volumePowerLossObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from time, dictionary and adjoint solver
        volumePowerLossObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~volumePowerLossObjectiveFunctionObject();


    // Member Functions

        //- read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate cost function value
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
