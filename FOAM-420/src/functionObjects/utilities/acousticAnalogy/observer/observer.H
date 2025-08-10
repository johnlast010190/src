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

Class
    Foam::observer

Description

    Stores observed pressure fluctuations


SourceFiles
    observer.C

\*---------------------------------------------------------------------------*/

#ifndef observer_H
#define observer_H

#include "db/IOstreams/Fstreams/OFstream.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tuple2/Tuple2.H"

#include "fvMesh/fvMesh.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "fields/Fields/scalarField/scalarField.H"

#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class observer Declaration
\*---------------------------------------------------------------------------*/

class observer
{
    const fvMesh& mesh_;

    //- Name of the sound observer
    word name_;

    //- Positions of the sound observer
    vectorField positions_;

    //- Current pressure fluctuations [Pa]
    scalarField pPrime_;
	scalarField pThickness_;
	scalarField pLoad_;
	scalarField pNonlinear_;
public:

    //- Runtime type information
    TypeName("observer");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            observer,
            word,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected observer
        static autoPtr<observer> New(
                                        const word& name,
                                        const fvMesh& mesh,
                                        const dictionary& dict
                                    );

    // Constructors

        // Construct from components
        observer
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Destructors
        virtual ~observer() {};


    // Member functions

        //- Return mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Set position of observer
        void setPositions(const vectorField& obsPos)
        {
            positions_ = obsPos;
        }

        //- Return position of observer
        virtual const vectorField& positions() const
        {
            return positions_;
        }

        //- Return position of observer
        label size() const
        {
            return positions().size();
        }

        //- Return fluctuation pressure
        const scalarField& pPrime() const
        {
            return pPrime_;
        }

		const scalarField& pThickness() const
        {
            return pThickness_;
        }

        scalarField& pThickness()
        {
            return pThickness_;
        }

        const scalarField& pLoad() const
        {
            return pLoad_;
        }

        scalarField& pLoad()
        {
            return pLoad_;
        }

        const scalarField& pNonlinear() const
        {
            return pNonlinear_;
        }

        scalarField& pNonlinear()
        {
            return pNonlinear_;
        }


        //- Return reference to fluctuation pressure
        scalarField& pPrime()
        {
            return pPrime_;
        }

        //- Set fluctuating pressure
        void pPrime(const scalarField& pPrime)
        {
            pPrime_ = pPrime;
        }

		void pThickness(const scalarField& pthickness)
        {
            pThickness_ = pthickness;
        }

        void pLoad(const scalarField& pload)
        {
            pLoad_ = pload;
        }

        void pNonlinear(const scalarField& pnonlinear)
        {
            pNonlinear_ = pnonlinear;
        }

        //- Return true to write on FO file
        virtual bool write() const
        {
            return true;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
