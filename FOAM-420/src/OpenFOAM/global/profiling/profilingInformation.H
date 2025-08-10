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
    (c) 2009-2016 Bernhard Gschaider
    (c) 2016-2107 OpenCFD Ltd.

Class
    Foam::profilingInformation

Description
    Code profiling information in terms of time spent, number of calls etc.

SourceFiles
    profilingInformation.C

\*---------------------------------------------------------------------------*/

#ifndef profilingInformation_H
#define profilingInformation_H

#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/strings/string/string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class profilingInformation;
class Ostream;

// Forward declaration of friend functions and operators
Ostream& operator<<(Ostream& os, const profilingInformation& info);


/*---------------------------------------------------------------------------*\
                    Class profilingInformation Declaration
\*---------------------------------------------------------------------------*/

class profilingInformation
{
    // Private Static Data Members

        //- Counter to generate the ids
        static label nextId_;

        //- Get a new ID and update the counter
        static label getNextId();

        //- Raise the next possible ID (to avoid ID-clashes during reading)
        static void raiseId(label maxVal);


    // Private Data Members

        //- Unique id to identify it
        const label id_;

        //- What this timer does
        const string description_;

        //- Pointer to the parent object (or self for top-level)
        profilingInformation* parent_;

        //- Nr of times this was called
        long calls_;

        //- Total time spent
        scalar totalTime_;

        //- Time spent in children
        scalar childTime_;

        //- Max memory usage on call.
        //  Only valid when the calling profiling has memInfo active.
        mutable int maxMem_;

        //- Is this information currently on the stack?
        mutable bool onStack_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        profilingInformation(const profilingInformation&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const profilingInformation&) = delete;


protected:

    // Friendship

        friend class profiling;


    // Constructors

        //- Construct null - only the master-element
        profilingInformation();


    // Member Functions

        //- Mark as being on the stack
        void push() const;

        //- Mark as being off the stack
        void pop() const;


        //- Write the profiling times, optionally with additional values
        //  Use dictionary format.
        Ostream& write
        (
            Ostream& os,
            const bool offset = false,
            const scalar elapsedTime = 0,
            const scalar childTime = 0
        ) const;


public:

    // Constructors

        //- Construct from components
        profilingInformation
        (
            const string& descr,
            profilingInformation* parent
        );


    //- Destructor
    ~profilingInformation();


    // Member Functions

      // Access

        inline label id() const
        {
            return id_;
        }


        inline const string& description() const
        {
            return description_;
        }


        inline profilingInformation& parent() const
        {
            return *parent_;
        }


        inline label calls() const
        {
            return calls_;
        }


        inline scalar totalTime() const
        {
            return totalTime_;
        }


        inline scalar childTime() const
        {
            return childTime_;
        }


        inline int maxMem() const
        {
            return maxMem_;
        }


        inline bool onStack() const
        {
            return onStack_;
        }


      // Edit

        //- Update it with a new timing information
        void update(const scalar elapsedTime);


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const profilingInformation& info
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
