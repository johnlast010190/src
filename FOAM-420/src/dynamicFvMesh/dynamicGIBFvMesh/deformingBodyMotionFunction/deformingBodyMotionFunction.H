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

Namespace
    Foam::deformingBodyMotionFunctions

Description
    Namespace for deforming-body motions


Class
    Foam::deformingBodyMotionFunction

Description
    Base class for defining interface motion for GIB

SourceFiles
    deformingBodyMotionFunction.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef deformingBodyMotionFunction_H
#define deformingBodyMotionFunction_H

#include "db/Time/Time.H"
#include "db/dictionary/dictionary.H"
#include "primitives/septernion/septernion.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class deformingBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class deformingBodyMotionFunction
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary DBMFCoeffs_;
        const Time& time_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        deformingBodyMotionFunction(const deformingBodyMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const deformingBodyMotionFunction&);


public:

    //- Runtime type information
    TypeName("deformingBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            deformingBodyMotionFunction,
            dictionary,
            (const fvMesh& mesh, const dictionary& SBMFCoeffs, const Time& runTime),
            (mesh, SBMFCoeffs, runTime)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        deformingBodyMotionFunction
        (
            const fvMesh&,
            const dictionary&,
            const Time&
        );


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<deformingBodyMotionFunction> New
        (
            const fvMesh&,
            const dictionary&,
            const Time&
        );


    //- Destructor
    virtual ~deformingBodyMotionFunction();


    // Member Functions

        //- Return the solid-body motion transformation septernion
//        virtual septernion transformation() const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;

        //- update
        virtual void update() = 0;

        //- Return the boundary velocity
        virtual tmp<vectorField> boundaryVelocity(const label&) const;

        //- Return the interface velocity
        virtual tmp<vectorField> interfaceVelocity(const label&) const;

        //- Return the interface velocity
        virtual tmp<vectorField> interfacePointsVelocity(const label&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
