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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::blockVertex

Description
    Define a block vertex.

SourceFiles
    blockVertex.C

\*---------------------------------------------------------------------------*/

#ifndef blockVertex_H
#define blockVertex_H

#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class blockVertex Declaration
\*---------------------------------------------------------------------------*/

class blockVertex
{

public:

    //- Runtime type information
    TypeName("blockVertex");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockVertex,
            Istream,
            (
                const dictionary& dict,
                const label index,
                const searchableSurfaces& geometry,
                Istream& is
            ),
            (dict, index, geometry, is)
        );


    // Constructors

        //- Construct null
        blockVertex();

        //- Clone function
        virtual autoPtr<blockVertex> clone() const;

        //- New function which constructs and returns pointer to a blockVertex
        static autoPtr<blockVertex> New
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            Istream&
        );

        //- Class used for the read-construction of
        //  PtrLists of blockVertex
        class iNew
        {
            const dictionary& dict_;
            const searchableSurfaces& geometry_;
            mutable label index_;

        public:

            iNew(const dictionary& dict, const searchableSurfaces& geometry)
            :
                dict_(dict),
                geometry_(geometry),
                index_(0)
            {}

            autoPtr<blockVertex> operator()(Istream& is) const
            {
                return blockVertex::New(dict_, index_++, geometry_, is);
            }
        };


    //- Destructor
    virtual ~blockVertex()
    {}


    // Member Functions

        virtual operator point() const = 0;

        //- Read vertex index with optional name lookup
        static label read(Istream&, const dictionary&);

        //- Write vertex index with optional name backsubstitution
        static void write(Ostream&, const label, const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
