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
    Foam::cuboid

Group
    grpRigidBodyDynamicsBodies

Description
    Specialization of rigidBody to construct a cuboid given the mass and
    lengths of the sides.

SourceFiles
    cuboidI.H
    cuboid.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_cuboid_H
#define RBD_cuboid_H

#include "bodies/rigidBody/rigidBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

/*---------------------------------------------------------------------------*\
                         Class cuboid Declaration
\*---------------------------------------------------------------------------*/

class cuboid
:
    public rigidBody
{
    // Private member data

        //- Lengths of the sides
        vector L_;


    // Private member functions

        // Calculate and return the inertia tensor
        inline symmTensor I(const scalar m, const vector& L) const;


public:

    //- Runtime type information
    TypeName("cuboid");


    // Constructors

        //- Construct from name, mass, centre of mass and lengths of the sides.
        inline cuboid
        (
            const word& name,
            const scalar m,
            const vector& c,
            const vector& L
        );

        //- Construct from dictionary
        inline cuboid
        (
            const word& name,
            const dictionary& dict
        );

        //- Return clone of this cuboid
        virtual autoPtr<rigidBody> clone() const;


    //- Destructor
    virtual ~cuboid();


    // Member Functions

        //- Return the lengths of the sides of the cuboid
        inline const vector& L() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodies/cuboid/cuboidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
