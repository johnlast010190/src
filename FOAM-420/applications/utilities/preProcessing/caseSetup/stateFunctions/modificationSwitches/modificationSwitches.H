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
    Foam::modificationSwitches

Description
    Simple container and access class for state modification switches

    resetInternalFields (true) overwrite internal fields
    resetBoundaryFields (true) overwrite boundary fields
                               (modified boundaries are always overwritten)
    resetSystemDicts (true) overwrite dictionaries in system folder
    resetConstDicts (true) overwrite dictionaries in constant folder
    resetBoundaryMesh (true) change boundary mesh
    deleteUnusedFields (true) delete fields that are not set by current state
    strictPatchNameChecking (true) exit with fatal error if search term
                                    does not find a match


SourceFiles
    modificationSwitchesI.H
    modificationSwitches.C

\*---------------------------------------------------------------------------*/

#ifndef modificationSwitches_H
#define modificationSwitches_H

#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class modificationSwitches Declaration
\*---------------------------------------------------------------------------*/

//- modification switch container
class modificationSwitches
{
    // Private data

        //- field and dictionary modification switches

            //- internal field reset switch
            Switch resetInternalFields_;

            //- boundary field reset switch
            Switch resetBoundaryFields_;

            //- system dictionary reset switch
            Switch resetSystemDicts_;

            //- const dictionary reset switch
            Switch resetConstDicts_;

            //- uniform dictionary reset switch
            Switch resetUniformDicts_;

            //- boundary mesh reset switch
            Switch resetBoundaryMesh_;

            //- switch for deleting fields not in current setup specification
            Switch deleteUnusedFields_;

            //- strict patch name checking switch
            Switch strictPatchNameChecking_;

            //- reuse existing dictionaries (defaults->existing->user)
            Switch reuseExistingDictionaries_;

            //- print boundary conditions to screen
            Switch printBoundaries_;

            //- Delete fields after each instanse initilisation (for multi-Instance only)
            Switch cleanInstanceFields_;


public:

    // Constructors

        //- Null constructor
        modificationSwitches();

        //- Construct from a dictionary,
        modificationSwitches
        (
            const dictionary& input
        );

        //- Construct from dictionary using provided switches as default
        modificationSwitches
        (
            const dictionary& input,
            const modificationSwitches& sw
        );

        //- Copy construct
        modificationSwitches
        (
            const modificationSwitches&
        );



    //- Destructor
    ~modificationSwitches();


    // Member Functions

        // Access

            //- const access to internal field reset switch
            inline bool resetInternalFields() const
            {
                return resetInternalFields_;
            }

            //- const access to boundary field reset switch
            inline bool resetBoundaryFields() const
            {
                return resetBoundaryFields_;
            }

            //- const access to boundary field reset switch
            inline const Switch& resetBoundaryFieldsRef() const
            {
                return resetBoundaryFields_;
            }

            //- const access to system dict modification switch
            inline bool resetSystemDicts() const
            {
                return resetSystemDicts_;
            }

            //- const access to constant dict modification switch
            inline bool resetConstDicts() const
            {
                return resetConstDicts_;
            }

            //- const access to constant dict modification switch
            inline bool resetUniformDicts() const
            {
                return resetUniformDicts_;
            }

            //- const access to boundary mesh modification switch
            inline bool resetBoundaryMesh() const
            {
                return resetBoundaryMesh_;
            }

            //- const access to switch for deleting unused fields
            inline bool deleteUnusedFields() const
            {
                return deleteUnusedFields_;
            }

            //- const access to strict patch name checking switch
            inline bool strictPatchNameChecking() const
            {
                return strictPatchNameChecking_;
            }

            //- const access to dictionary reuse switch
            inline bool reuseExistingDictionaries() const
            {
                return reuseExistingDictionaries_;
            }

            //- const access to dictionary reuse switch
            inline bool printBoundaries() const
            {
                return printBoundaries_;
            }

            inline bool cleanInstanceFields() const
            {
                return cleanInstanceFields_;
            }


        // Check

        // Edit

            //- merge dictionary with existing switches
            void merge(const dictionary&);

        // Write


    // Member Operators

        void operator=(const modificationSwitches&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
