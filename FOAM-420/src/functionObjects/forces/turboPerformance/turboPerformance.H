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
    (c) D. Boger, B. Lewis M. Auvinen, H. Nilsson

Class
    Foam::functionObjects::turboPerformance

Description


Authors
    Mikko Auvinen

Modifications
    Ported to OpenFoam-1.6-ext using patch by David Bogar - Penn State, ARL

    Added capability to analyse both pumps and turbines using a flag in the
    controlDict by Bryan Lewis - Penn State

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    turboPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turboPerformance_H
#define functionObjects_turboPerformance_H

#include "turboPerformance/fluidPower.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class turboPerformance Declaration
\*---------------------------------------------------------------------------*/

class turboPerformance
:
    public fluidPower
{
    // Private data

       // - Rotational velocity (in rad/s)
        vector omega_;

        //- moving wall sections
        const label movingWallSectors_;

    // Private member functions

        //- Disallow default bitwise copy construct
        turboPerformance(const turboPerformance&);

        //- Disallow default bitwise assignment
        void operator=(const turboPerformance&);


protected:

    //- Output file header information
    virtual void writeFileHeader(Ostream& os) const;

public:

    //- Runtime type information
    TypeName("turboPerformance");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        turboPerformance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~turboPerformance();


    // Member Functions


        //- Execute, currently does nothing
        virtual bool execute();

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Calculate the forces and write
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
