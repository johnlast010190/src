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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::extrudeModel

Description
    Top level extrusion model class

SourceFiles
    extrudeModel.C

\*---------------------------------------------------------------------------*/

#ifndef extrudeModel_H
#define extrudeModel_H

#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/point.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class extrudeModel Declaration
\*---------------------------------------------------------------------------*/

class extrudeModel
{
protected:

    // Protected data

        label nLayers_;

        const scalar expansionRatio_;

        const dictionary& dict_;

        const dictionary& coeffDict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        extrudeModel(const extrudeModel&);

        //- Disallow default bitwise assignment
        void operator=(const extrudeModel&);


public:

    //- Runtime type information
    TypeName("extrudeModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            extrudeModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from dictionary
        extrudeModel(const word& modelType, const dictionary&);


    // Selectors

        //- Select null constructed
        static autoPtr<extrudeModel> New(const dictionary&);


    //- Destructor
    virtual ~extrudeModel();


    // Member Functions

        // Access

            label nLayers() const;

            scalar expansionRatio() const;


    // Member Operators

        //- Helper: calculate cumulative relative thickness for layer.
        //  (layer=0 -> 0; layer=nLayers -> 1)
        scalar sumThickness(const label layer) const;

        virtual point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
