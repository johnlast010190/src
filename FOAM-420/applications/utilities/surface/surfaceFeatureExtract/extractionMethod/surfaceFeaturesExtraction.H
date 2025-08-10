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
    (c) 2017 OpenCFD Ltd.

Namespace
    Foam::surfaceFeaturesExtraction

Description
    Namespace for run-time selectable surface feature extraction methods.

Class
    Foam::surfaceFeaturesExtraction::method

Description
    Abstract base for run-time selectable surface feature extraction methods.

SourceFiles
    surfaceFeaturesExtraction.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFeaturesExtraction_method_H
#define surfaceFeaturesExtraction_method_H

#include "triSurface/surfaceFeatures/surfaceFeatures.H"
#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{

/*---------------------------------------------------------------------------*\
              Class surfaceFeaturesExtraction::method Declaration
\*---------------------------------------------------------------------------*/

class method
{
protected:

    scalar includedAngle_;
    Switch geometricTestOnly_;

    //- Construct null
    method();


public:

    //- Runtime type information
    ClassNameNoDebug("method");


    // Constructors

        //- Construct from dictionary
        method(const dictionary& dict);


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        method,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<method> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~method();


    // Member Functions

        //- The included angle, if set
        inline scalar includedAngle() const
        {
            return includedAngle_;
        }

        //- Use geometric test only
        inline Switch geometricTestOnly() const
        {
            return geometricTestOnly_;
        }

        //- Extracted features
        virtual autoPtr<surfaceFeatures> features
        (
            const triSurface& surf
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFeaturesExtraction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
