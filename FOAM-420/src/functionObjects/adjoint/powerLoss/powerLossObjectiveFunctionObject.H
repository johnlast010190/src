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
    powerLossObjectiveFunctionObject

Description
    Calculates the value for power losses between inlet and outlet.

SourceFiles
    powerLossObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef powerLossObjectiveFunctionObject_H
#define powerLossObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class powerLossObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class powerLossObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected member data

        //- Total objective patch inlet area
        scalar areaObjPatchInlet_;

        //- Total objective patch outlet area
        scalar areaObjPatchOutlet_;

        //- Switch on adiabatic reporting for compressible flow (default false)
        Switch adiabatic_;

        //- Switch on tangential boundary stress energy transport
        //  Defaults to false
        Switch stressWork_;

        //- Formulation of objective for compressible cases
        word compressibleForm_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        powerLossObjectiveFunctionObject
        (
            const powerLossObjectiveFunctionObject&
        );
        void operator=(const powerLossObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("powerLoss");


    // Constructors

        //- Construct from time, dictionary, adjoint solver and use old-style
        //  adjoint file-format
        powerLossObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from time, dictionary and adjoint solver
        powerLossObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~powerLossObjectiveFunctionObject();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Execute the averaging
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
