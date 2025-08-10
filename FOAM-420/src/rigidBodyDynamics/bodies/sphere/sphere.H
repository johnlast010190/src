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
    Foam::sphere

Group
    grpRigidBodyDynamicsBodies

Description
    Specialization of rigidBody to construct a sphere given the mass and radius.

SourceFiles
    sphereI.H
    sphere.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_sphere_H
#define RBD_sphere_H

#include "bodies/rigidBody/rigidBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

/*---------------------------------------------------------------------------*\
                         Class sphere Declaration
\*---------------------------------------------------------------------------*/

class sphere
:
    public rigidBody
{
    // Private member data

        //- Radius
        scalar r_;


    // Private member functions

        // Calculate and return the inertia tensor
        inline symmTensor I(const scalar m, const scalar r) const;


public:

    //- Runtime type information
    TypeName("sphere");


    // Constructors

        //- Construct from name, mass, centre of mass and radius
        inline sphere
        (
            const word& name,
            const scalar m,
            const vector& c,
            const scalar r
        );

        //- Construct from dictionary
        inline sphere
        (
            const word& name,
            const dictionary& dict
        );

        //- Return clone of this sphere
        virtual autoPtr<rigidBody> clone() const;


    //- Destructor
    virtual ~sphere();


    // Member Functions

        //- Return the radius of the sphere
        inline scalar r() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodies/sphere/sphereI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
