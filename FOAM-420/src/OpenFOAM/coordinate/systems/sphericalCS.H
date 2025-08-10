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
    (c) 2010-2016 Esi Ltd.

Class
    Foam::sphericalCS

Description
    A spherical coordinate system (r-theta-phi).
    The coordinate system angles theta and phi are always in radians.

SourceFiles
    sphericalCS.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalCS_H
#define sphericalCS_H

#include "coordinate/systems/coordinateSystem.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{
/*---------------------------------------------------------------------------*\
                         Class spherical Declaration
\*---------------------------------------------------------------------------*/

class spherical
:
    public coordinateSystem
{
    // Private data members

protected:

    // Protected Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual vector localToGlobal(const vector&, bool translate) const;

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual vector globalToLocal(const vector&, bool translate) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField&,
            bool translate
        ) const;


public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct null
        spherical();

        //- Copy construct
        spherical(const spherical& csys) = default;

        //- Move construct
        spherical(spherical&& csys) = default;

        //- Copy construct from another coordinateSystem type
        explicit spherical(const coordinateSystem& csys);

        //- Move construct from another coordinateSystem type
        explicit spherical(coordinateSystem&& csys);

        //- Move construct from autoPtr of another coordinateSystem type
        explicit spherical(autoPtr<coordinateSystem>&& csys);

        //- Construct copy with a different name
        spherical
        (
            const word& name,
            const coordinateSystem&
        );

        //- Construct from origin and rotation
        spherical
        (
            const point& origin,
            const coordinateRotation&
        );

        //- Construct from origin and 2 axes
        spherical
        (
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from origin and 2 axes
        spherical
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary
        spherical(const word& name, const dictionary&);

        //- Construct from dictionary without a name
        explicit spherical(const dictionary&);

        //- Construct from dictionary with optional subDict lookup.
        //
        //  \param dictName If non-empty, the sub-dictionary to use.
        spherical(const dictionary& dict, const word& dictName);

        //- Return clone
        virtual autoPtr<coordinateSystem> clone() const
        {
            return autoPtr<coordinateSystem>::NewFrom<spherical>(*this);
        }

        //- Convert from Cartesian (to Spherical)
        static inline vector fromCartesian(const vector& v)
        {
            return
                vector
                (
                    sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z()),
                    atan2(hypot(v.x(), v.y()), v.z()),
                    atan2(v.y(),v.x())
                );
        }

        //- Convert to Cartesian (from Spherical)
        static inline vector toCartesian(const vector& v)
        {
            return
                vector
                (
                    v.x()*sin(v.y())*cos(v.z()),
                    v.x()*sin(v.y())*sin(v.z()),
                    v.x()*sin(v.y())
                );
        }


    //- Destructor
    virtual ~spherical();

    // Member Functions

        //- Treat the rotation tensor as non-uniform
        virtual bool uniform() const
        {
            return false;
        }

        //- Position-dependent rotation tensors at multiple points
        using coordinateSystem::R;

        //- Position-dependent rotation tensor at a single point
        //- \return tensor
        virtual tensor R(const point& global) const;

    // Member Operators

        //- Copy assignment
        spherical& operator=(const spherical&) = default;

        //- Move assignment
        spherical& operator=(spherical&&) = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSystem
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
