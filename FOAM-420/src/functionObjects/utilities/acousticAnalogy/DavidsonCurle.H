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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::functionObjects::DavidsonCurle

Group
    grpUtilitiesFunctionObjects

Description
    Implementation of the DavidsonCurle analogy for incompressible flow as
    described in:
    \verbatim
        Johan Larsson, Lars Davidson, Magnus Olsson, Lars-Erik Eriksson (2004).
        Aeroacoustic Investigation of an Open Cavity at Low Mach Number.
        AIAA Journal, Vol. 42, No. 12.
    \endverbatim

    Example of function object specification:
    \verbatim
    DavidsonCurle1
    {
        type                    DavidsonCurle;
        functionObjectLibs      ("libAcousticAnalogy.so");
        log                     true;

        patches                 ( cylinder1 cylinder2 );
        cellZone                cellZoneName;
        rhoRef                  1.225;
        cRef                    343;

        observers
        {
            microphone-A
            {
                position        (-0.4760595 1.58962725 0);
            }
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description | Required | Default value
        type         | type name: DavidsonCurle | yes |
        log          | write acoustic data to standard output  | no | no
        patches      | patches included in the acoustic calculation    | yes |
        cellZone     | name of cell zone for acoustic volume integrals | no |
        rhoRef       | reference density   | yes |
        cRef         | reference speed of sound  | yes |
        pName        | pressure field name | no  | p
        UName        | velocity field name | no  | U
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    DavidsonCurle.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DavidsonCurle_H
#define functionObjects_DavidsonCurle_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "acousticAnalogy/soundObserver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class DavidsonCurle Declaration
\*---------------------------------------------------------------------------*/

class DavidsonCurle
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Initialised flag
        bool initialised_;

        //- Patches to sample
        labelHashSet patches_;

        //- Name of cellZone
        word cellZoneName_;

        //- ID of cellZone
        label cellZoneID_;

        //- Name of pressure field (optional)
        word pName_;

        //- Name of velocity field (optional)
        word UName_;

        //- Reference density needed for incompressible cases
        scalar rhoRef_;

        //- Reference speed of sound
        scalar cRef_;

        //- List of observers
        List<SoundObserver> observers_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Helper function to write acoustic data
        void writeDavidsonCurle();

        //- Initialise the fields
        void initialise();

        //- Return pressure
        tmp<volScalarField> p() const;

        //- Return pressure time derivative
        tmp<volScalarField> dpdt() const;

        //- Return Lighthill tensor
        tmp<volTensorField> Tij() const;

        //- Return first time derivative of Lighthill tensor
        tmp<volTensorField> dTijdt() const;

        //- Return second time derivative of Lighthill tensor
        tmp<volTensorField> d2Tijdt2() const;

        //- Disallow default bitwise copy construct
        DavidsonCurle(const DavidsonCurle&);

        //- Disallow default bitwise assignment
        void operator=(const DavidsonCurle&);


public:

    //- Runtime type information
    TypeName("DavidsonCurle");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        DavidsonCurle
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~DavidsonCurle();


    // Member Functions

        //- Calculate the acoustic data
        virtual void calculate();

        //- Read the DavidsonCurle data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the acoustic data
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
