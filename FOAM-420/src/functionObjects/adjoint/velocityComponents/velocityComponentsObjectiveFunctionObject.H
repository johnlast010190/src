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
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::objectiveFunctionObject::velocityComponentsFunctionObject

Description


SourceFiles
    velocityComponentsObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef velocityComponentsObjectiveFunctionObject_H
#define velocityComponentsObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

class velocityComponentsObjectiveFunctionObject: public objectiveFunctionObject
{

protected:

    vector axis_;

    vector origin_;

    vector swirlDirection_;

    bool volumeBased_;

    scalar objectiveVolume_;
    scalar objectiveArea_;

    //-Actual values
    class velocityData{
    public:
        scalar target_;
        scalar weights_;
        scalar deviation_;
        bool active_;

        velocityData()
        {
            target_ = 0.;
            weights_ = 0.;
            deviation_ = -1;
            active_ = false;
        }
    };

    velocityData axialVel_, peripheralVel_, radialVel_, peripheralOverAxialVel_;

    //- Write the header of the output file in postProcessing
    virtual void writeFileHeader(Ostream&);

    bool isAnActivePatch(const label& pI);

    void setVariables(const dictionary& objectiveDict);
    void calculateObjectiveVolume();
    void calculateObjectiveArea();

    scalar getVolumeAxialDeviation();
    scalar getVolumePeripheralDeviation();
    scalar getVolumeRadialDeviation();
    scalar getVolumePeripheralOverAxialDeviation();

    scalar getSurfaceAxialDeviation();
    scalar getSurfacePeripheralDeviation();
    scalar getSurfaceRadialDeviation();
    scalar getSurfacePeripheralOverAxialDeviation();

private:

    //- Disallow default bitwise copy construct and assignment
    velocityComponentsObjectiveFunctionObject
    (
        const velocityComponentsObjectiveFunctionObject&
    );
    void operator=(const velocityComponentsObjectiveFunctionObject&);


public:

    TypeName("velocityComponents");

    //- Constructors
        velocityComponentsObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        velocityComponentsObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

    //- Destructor
    virtual ~velocityComponentsObjectiveFunctionObject();

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
