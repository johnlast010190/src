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
    turbulentNoiseObjectiveFunctionObject

Description
    Noise approximated by the square of nut integrated over a volume.

    "Noise reduction in car aerodynamics using a surrogate objective function
    and the continuous adjoint method with wall functions"
    E.M. Papoutsis-Kiachagias*,
    N.   Magoulas*,
    J.   Mueller**,
    C.   Othmer**,
    K.C. Giannakoglou*,
    *National Technical University of Athens, School of Mechanical Engineering,
     Parallel CFD & Optimization Unit,
     P.O. Box 64069, 15710 Athens, Greece
    ** Volkswagen AG, CAE Methods, Group Research,
       Letter Box 1777, D-38436 Wolfsburg, Germany

SourceFiles
    turbulentNoiseObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentNoiseObjectiveFunctionObject_H
#define turbulentNoiseObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class turbulentNoiseObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class turbulentNoiseObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        wordList  aeroAcousticsZoneNames_;
        labelList aeroAcousticsZoneLabels_;


private:

    // Private Member Functions

        labelList findZoneIDs();

        //- Disallow default bitwise copy construct and assignment
        turbulentNoiseObjectiveFunctionObject
        (
            const turbulentNoiseObjectiveFunctionObject&
        );
        void operator=(const turbulentNoiseObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("turbulentNoise");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        turbulentNoiseObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        turbulentNoiseObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~turbulentNoiseObjectiveFunctionObject();


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
