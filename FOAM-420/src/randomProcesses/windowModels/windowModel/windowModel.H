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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::windowModel

Description
    Base class for windowing models

SourceFiles
    noiseFFT.C

\*---------------------------------------------------------------------------*/

#ifndef windowModel_H
#define windowModel_H

#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/Fields/scalarField/scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class windowModel Declaration
\*---------------------------------------------------------------------------*/

class windowModel
:
    public scalarField
{

protected:

    // Protected Data

        //- Number of overlap samples per window
        label nOverlapSamples_;

        //- Number of windows
        label nWindow_;


public:

    //- Runtime type information
    TypeName("windowModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        windowModel,
        dictionary,
        (
            const dictionary& dict,
            const label nSamples
        ),
        (dict, nSamples)
    );


    //- Construct from dictionary
    windowModel(const dictionary& dict, const label nSamples);


    // Selectors

        //- Return a reference to the selected window model
        static autoPtr<windowModel> New
        (
            const dictionary& dict,
            const label nSamples
        );


    //- Destuctor
    virtual ~windowModel();


    // Public Member Functions

        //- Return the number of samples in the window
        label nSamples() const;

        //- Return the number of windows
        label nWindow() const;

        //- Return the total number of windows for a given number of samples
        label nWindowsTotal(label nSamplesTotal) const;

        //- Validate that the window is applicable to the data set size, and
        //  return the number of required data points
        label validate(label n);

        //- Return the windowed data
        template<class Type>
        tmp<Field<Type>> apply
        (
            const Field<Type>& fld,
            const label windowI
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "windowModels/windowModel/windowModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
