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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::lduMesh

Description
    Abstract base class for meshes which provide LDU addressing for the
    construction of lduMatrix and LDU-solvers.

\*---------------------------------------------------------------------------*/

#ifndef lduMesh_H
#define lduMesh_H

#include "matrices/lduMatrix/lduAddressing/lduAddressing.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterfacePtrsList.H"
#include "db/typeInfo/typeInfo.H"
#include "db/IOstreams/IOstreams/InfoProxy.H"
#include "primitives/Vector/lists/vectorList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;


// Forward declaration of friend functions and operators

class lduMesh;

Ostream& operator<<(Ostream&, const InfoProxy<lduMesh>&);


/*---------------------------------------------------------------------------*\
                           Class lduMesh Declaration
\*---------------------------------------------------------------------------*/

class lduMesh
{

        //- Identifier of boundary cells
        boolList isBoundaryCell_ = boolList(0);

public:

    //- Runtime type information
    TypeName("lduMesh");


    // Constructors

    //- Destructor
    virtual ~lduMesh()
    {}


    // Member Functions

        // Access

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const = 0;

            //- Return the object registry
            virtual const objectRegistry& thisDb() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const = 0;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const = 0;

            //- Return communicator used for parallel communication
            virtual label comm() const = 0;

            //- Helper: reduce with current communicator
            template<class T, class BinaryOp>
            void reduce
            (
                T& Value,
                const BinaryOp& bop
            ) const;

            //- isBoundaryCell
            virtual const boolList& isBoundaryCell() const
            {
                return isBoundaryCell_;
            }

            //- isBoundaryCell
            virtual boolList& isBoundaryCell()
            {
                return isBoundaryCell_;
            }

        // Info

            //- Return info proxy.
            //  Used to print mesh information to a stream
            InfoProxy<lduMesh> info() const
            {
                return *this;
            }


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const InfoProxy<lduMesh>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshes/lduMesh/lduMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
