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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::windowModels::Hanning

Description
    Hanning window

    The window is described by the function
    \f[
        wf = (1 - \alpha) - \alpha cos(2 \pi t/m);
    \f]

    Where:
    \vartable
        \alpha  | Coefficient with a default value of 0.5
        t       | time
        m       | window width
    \endvartable

    The window can be further manipulated by the controls:
    - \c symmetric: force the window to be symmetric
    - \c extended: extend the window by 1 element at start and end to produce
      non-zero values at the start and end positions.  Note: window is
      normalised to preserve energy content

SourceFiles
    Hanning.C

\*---------------------------------------------------------------------------*/

#ifndef Hanning_H
#define Hanning_H

#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "windowModels/windowModel/windowModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace windowModels
{

/*---------------------------------------------------------------------------*\
                           Class Hanning Declaration
\*---------------------------------------------------------------------------*/

class Hanning
:
    public windowModel
{

protected:

    // Protected Member Data

        //- Symmetric switch
        bool symmetric_;

        //- Extended switch
        bool extended_;

        //- Window coefficient, default = 0.5
        scalar alpha_;


public:

    //- Runtime type information
    TypeName("Hanning");


    //- Construct from dictionary
    Hanning(const dictionary& dict, const label nSamples);

    //- Destuctor
    virtual ~Hanning();


    // Public Member Functions

        //- Return the symmetric flag
        bool symmetric() const;

        //- Return the extended flag
        bool extended() const;

        //- Return the window coefficient
        scalar alpha() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace windowModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
