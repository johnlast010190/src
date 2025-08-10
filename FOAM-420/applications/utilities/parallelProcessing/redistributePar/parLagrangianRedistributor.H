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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::parLagrangianRedistributor

Description
    Lagrangian field redistributor.

    Runs in parallel. Redistributes from fromMesh to toMesh.

SourceFiles
    parLagrangianRedistributor.C
    parLagrangianRedistributorRedistributeFields.C

\*---------------------------------------------------------------------------*/

#ifndef parLagrangianRedistributor_H
#define parLagrangianRedistributor_H

#include "containers/Lists/PtrList/PtrList.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistributePolyMesh;
class mapDistributeBase;
class IOobjectList;
class passiveParticleCloud;

/*---------------------------------------------------------------------------*\
                    Class parLagrangianRedistributor Declaration
\*---------------------------------------------------------------------------*/

class parLagrangianRedistributor
{
    // Private data

        //- Source mesh reference
        const fvMesh& srcMesh_;

        //- Destination mesh reference
        const fvMesh& tgtMesh_;

        //- Distribution map reference
        const mapDistributePolyMesh& distMap_;

        //- For every src cell the target processor
        labelList destinationProcID_;

        //- For every src cell the target cell
        labelList destinationCell_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        parLagrangianRedistributor(const parLagrangianRedistributor&);

        //- Disallow default bitwise assignment
        void operator=(const parLagrangianRedistributor&);


public:

    // Constructors

        //- Construct from components
        parLagrangianRedistributor
        (
            const fvMesh& srcMesh,
            const fvMesh& tgtMesh,
            const label nOldCells,
            const mapDistributePolyMesh& distMap
        );


    // Member Functions

        //- Find all clouds (on all processors) and for each cloud all
        //  the objects. Result will be synchronised on all processors
        static void findClouds
        (
            const fvMesh&,
            wordList& cloudNames,
            List<wordList>& objectNames
        );

        //- Redistribute and write lagrangian positions
        autoPtr<mapDistributeBase> redistributeLagrangianPositions
        (
            passiveParticleCloud& cloud
        ) const;

        //- Read, redistribute and write lagrangian positions
        autoPtr<mapDistributeBase> redistributeLagrangianPositions
        (
            const word& cloudName
        ) const;

        //- Pick up any fields of a given type
        template<class Type>
        static wordList filterObjects
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Read, redistribute and write all/selected lagrangian fields
        template<class Type>
        void redistributeLagrangianFields
        (
            const mapDistributeBase& map,
            const word& cloudName,
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        ) const;

        //- Read, redistribute and write all/selected lagrangian fieldFields
        template<class Type>
        void redistributeLagrangianFieldFields
        (
            const mapDistributeBase& map,
            const word& cloudName,
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        ) const;

        //- Read and store all fields of a cloud
        template<class Container>
        static void readLagrangianFields
        (
            const passiveParticleCloud& cloud,
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Redistribute and write stored lagrangian fields
        template<class Container>
        void redistributeStoredLagrangianFields
        (
            const mapDistributeBase& map,
            passiveParticleCloud& cloud
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "parLagrangianRedistributorRedistributeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
