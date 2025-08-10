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
    (c) 2020-2020 Esi Ltd.

Class
    Foam::meshStatistics

Description
    Storage of mesh quality field and frequency statistics

SourceFiles
    meshStatistics.C

\*---------------------------------------------------------------------------*/

#ifndef meshStatistics_H
#define meshStatistics_H

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class meshStatistics Declaration
\*---------------------------------------------------------------------------*/

class meshStatistics
{
        //- Quality field name
        const word name_;

        //- Mesh statistics reporting dictionary
        const dictionary dict_;

        //- Quality volume field.
        volScalarField quality_;

        //- Quality frequency data
        labelList frequency_;

        //- Number of data samples
        label nSamples_;

        //- Ranges for frequency sampling
        scalar minRange_;
        scalar binSize_;

        //- Min mesh quality metric
        scalar minVal_;

        //- Max mesh quality metric
        scalar maxVal_;

        //- Average mesh quality metric
        scalar aveVal_;

        //- Standard deviation
        scalar stdDev_;

        //- Whether histograms have been calculated
        bool active_;

public:

    ClassName("meshStatistics");

    // Constructors

        //- Construct from mesh
        meshStatistics
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& dict
        );


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        const volScalarField& field() const
        {
            return quality_;
        }

        volScalarField& field()
        {
            return quality_;
        }

        void calcStats
        (
            const List<scalar>& metrics
        );

        void writeStats();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
