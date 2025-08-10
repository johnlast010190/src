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
    (c) 2010-2011, Esi Ltd

Class
    Foam::functionObjects::zoneForces

Group
    grpForcesFunctionObjects

Description
    Reporting of zones.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    zoneForces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_zoneForces_H
#define functionObjects_zoneForces_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class zoneForces Declaration
\*---------------------------------------------------------------------------*/

class zoneForces
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        label nAveragingSteps_;

        bool verbose_;

        bool logToFile_;

        //- Sorted list of porous zones
        List<word> porousZoneNames_;

        //- mass flux file ptr
        autoPtr<OFstream > massFluxFilePtr_;

        //- pressure drop file ptr
        autoPtr<OFstream > pressureDropFilePtr_;

        //- Result values
        List<labelList > boundaryFaces_;
        List<scalar > massFlux_;
        List<vector > force_;
        List<scalar > pressureDrop_;

        label averagingIndex_;

    // Private Member Functions

        //- return sum of mag of face flux
        scalar zoneMassFlux(const labelList& faces);

        //- return sum of faces areas projected onto force direction
        scalar zoneNormalForce
        (
            const labelList& faces,
            const vector& forces
        );

        //- Disallow default bitwise copy construct
        zoneForces(const zoneForces&);

        //- Disallow default bitwise assignment
        void operator=(const zoneForces&);

public:

    //- Runtime type information
    TypeName("zoneForces");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        zoneForces
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~zoneForces();


    // Member Functions

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if rhoName is specified otherwise rhoRef
        tmp<volScalarField> getRho() const;

        //- update values
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the zoneForces
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
