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
    Foam::sensorTypes::set

Description
    Sensor type based on mean/min/max value of vol<Type>Field in a given cell
    set or cell zone.

    \verbatim
        setType    <type>; //available: cellSet/cellZone;
        <type>     <setName>; // i.e. cellSet   c0; cellZone  z0;
        sensorOp   <operation>; //available: min/max/mean
    \endverbatim

SourceFiles
    set.C

\*---------------------------------------------------------------------------*/

#ifndef set_H
#define set_H

#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{

/*---------------------------------------------------------------------------*\
                           Class set Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class set
:
    public sensor<Type>
{
    // Private data

        //- set type
        //  currently supported: cellSet, cellZone
        const word setType_;

        //- sensor type
        const word sensorOp_;


public:

    // Runtime type information
    TypeName("set");


    // Constructors

        //- Construct from components
        set
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy constructor
        set(const set<Type>& cv);


    //- Destructor
    virtual ~set();


    // Member Functions

        //- Return set value
        Type value(const GeometricField<Type, fvPatchField, volMesh>& volField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sensorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/sensor/set/set.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
