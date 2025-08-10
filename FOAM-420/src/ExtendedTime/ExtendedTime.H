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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::ExtendedTime

Description
    Class to control time during OpenFOAM simulations that is also the
    top-level objectRegistry.

SourceFiles
    ExtendedTime.C

\*---------------------------------------------------------------------------*/


#ifndef ExtendedTime_H
#define ExtendedTime_H

#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ExtendedTime Declaration
\*---------------------------------------------------------------------------*/

class ExtendedTime
:
    public Time
{

private:

    // Private data

        //- The start time of the sensitivities computation interval
        scalar initialTime_;

        //- The end time of the sensitivities computation interval
        scalar finalTime_;

        wordList functionObjectsInReverse_;

        label reverseTimestep_;

        scalar primalDeltaT_;

        scalar reverseDeltaT_;

        scalar reverseDeltaTSave_;

        // Checkpoints

            //- Lookup table of checkpointing fields
            wordList cpTable_;

            //- Number of checkpoints
            label nCP_;

            //- Store checkpoints at the end of the primal
            Switch storeCPs_;

            //- List of pointers to checkpoints
            PtrList<objectRegistry> checkpoint_;

            //- Time value of the checkpoints
            List<scalar>  cpTimeValue_;

            //- deltaTvalue of the checkpoints
            List<scalar>  cpDeltaTValue_;

            //- deltaT0value of the checkpoints
            List<scalar>  cpDeltaT0Value_;

            //- Time index of the checkpoints
            List<label>  cpTimeIndex_;

            //- Level of the checkpoint
            List<label>  cpLevel_;

            //- Pointer to the checkpoints in the database
            List<label>  cpAddressing_;


    // Private member functions

        // Time controls

            //- Return true if run should continue,
            bool inverseRun() const;

            //- Return true if run should continue and if so increment time
            bool inverseLoop();

            //- Return true if end of run,
            bool inverseEnd() const;

        // Checkpointing controls

            //- Initialize checkpoints
            void init();

            //- Add a checkpoint at current time
            void addCP();

            //- Store fields in the specified checkpoint
            void storeCP(label);

            //- Check whether the checkpoint of the current time is available
            bool cpIsAvailable();

            //- Restore the checkpoint of the current time
            void restoreCurrentTimeCP();

            //- Restore the checkpoint from the latest available in memory
            void restoreLatestTimeCP();

            //- Restore the specified checkpoint
            void restoreCP(label);

            //- Copy current geometric fields to checkpoint
            template
            <
                class Type,
                template<class> class PatchField,
                class GeoMesh
            >
            void addCurrentFieldsToCP(label);

            //- Copy current geometric fields to checkpoint
            template
            <
                class Type,
                template<class> class PatchField,
                class GeoMesh
            >
            void restoreFieldsFromCP(label);

            bool isValidTimeStepInReverse(const scalar);

public:

    TypeName("ExtendedTime");

    // Constructors

        //- Construct given name of dictionary to read and argument list
        ExtendedTime
        (
            const word& name,
            const argList& args,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given name of dictionary to read, rootPath and casePath
        ExtendedTime
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given dictionary, rootPath and casePath
        ExtendedTime
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given rootPath and casePath
        ExtendedTime
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

    //- Destructor
    virtual ~ExtendedTime();


    // Member functions

        //- Return true if run should continue and if so decrease time
        bool loop();

        //- Return true if run should continue and if so decrease time
        bool backloop();

        //-
        bool forloop();


    // Member operators

        //- Set deltaT to that specified and increment time via operator++()
        ExtendedTime& operator-=(const dimensionedScalar&);

        //- Set deltaT to that specified and increment time via operator++()
        ExtendedTime& operator-=(const scalar);

        ExtendedTime& operator--();

        ExtendedTime& operator--(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ExtendedTimeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
