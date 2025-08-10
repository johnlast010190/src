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
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::distributeField

SourceFiles
    distributeField.C

Authors
    Daniel Rettenmaier   <rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising       <deising@mma.tu-darmstadt.de>
    Christian Kunkelmann (formerly TTD, TU Darmstadt)
    All rights reserved.

Description
    Distributes a field starting with a configuration marked by a marker set
    to a destination set by averaging field information to neighbor cells.
    Additionally the distance of the current cells to the start set can be
    calculated.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier (main developer).

    Method Development and Intellectual Property :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising     <deising@mma.tu-darmstadt.de>
        Holger Marschall   <marschall@mma.tu-darmstadt.de>
        Dieter Bothe       <bothe@mma.tu-darmstadt.de>
        Cameron Tropea     <ctropea@sla.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef distributeField_H
#define distributeField_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class distributeField Declaration
\*---------------------------------------------------------------------------*/

class distributeField
{

public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        distributeField(){};

    //- Destructor
        virtual ~distributeField(){};

    // Member Functions

        // distributes a volume field starting with a current distribution. The field
        // information will be transported to neighbor cells where the data is averaged
        // with the information of other neighbors. Several loops may be performed until
        // the information is transported to all destination cells. Optionally the
        // field can be normalized if it is a normal vector
        template <typename T>
        void distributeVolField
        (
            GeometricField<T, fvPatchField, volMesh>& field,   //field to distribute
            volScalarField& currField,         //numbered field with current distribution status
            const volScalarField& destField,   //numbered field with width of the destination field
            int maxLoops,
            int minLoops,
            bool normalizeField
        ) const;


        // distributes a normal volume vector field and calculates the distance of the
        // current cells to the start-configuration. Hereby the cells of the start-configuration
        // already contain the normal and the distance information.
        // The information will be transported to neighbor cells where the data is averaged
        // with the information of other neighbors. Several loops may be performed until
        // the information is transported to all destination cells. Optionally the
        // field can be normalized if it is a normal vector
        //
        // distance to interface
        //  > 0: liquid side
        //  < 0: vapor side
        void calculateDistance
        (
            boolList& currField,        //list numbered as field holding distance
            volScalarField& distance,   //information distance to interface
            volVectorField& nHatv,      //nHat must be propagated with distance
            const boolList& destField,  //minimum field to cover
            scalar maxDistance,
            int maxLoops,
            int minLoops
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interface/distributeField/distributeFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
