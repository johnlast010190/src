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
    (c) 2017 Esi Ltd.

Class
    Foam::sensorTypes::timeValue

Description
    Sensor type based on simulation time.


SourceFiles
    timeValue.C

\*---------------------------------------------------------------------------*/

#ifndef timeValue_H
#define timeValue_H

#include "cfdTools/general/sensor/sensor/sensor.H"
#include "sets/topoSetSource/topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{

/*---------------------------------------------------------------------------*\
                           Class timeValue Declaration
\*---------------------------------------------------------------------------*/

class timeValue
:
    public sensor<scalar>
{
    // Private data

        //- patch name
        word patchName_;


public:

    // Runtime type information
    TypeName("timeValue");


    // Constructors

        //- Construct from components
        timeValue
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy constructor
        timeValue(const timeValue& cv);


    //- Destructor
    virtual ~timeValue();


    // Member Functions

        //- Return timeValue value
        virtual tmp<Field<scalar>> valueField() const;

        //- Return timeValue value
        virtual scalar value() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sensorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
