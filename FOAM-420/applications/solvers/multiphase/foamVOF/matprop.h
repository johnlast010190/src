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
    (c) ESI Ltd, 2019.

 Description :

  This file declares the objects and functions for material properties used
  in FoamVOF model implementation. FoamVOF is a general-purpose multi-fluid
  module developed under the framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/

#ifndef __MATPROP_H
#define __MATPROP_H

namespace Foam {

class MatProp {
    protected:
        scalar cp1_;
        scalar cp2_;
        scalar k1_;
        scalar k2_;
        scalar rho1_;
        scalar rho2_;
        scalar Sch_;
        scalar Pr_;
        scalar psat_;
    public:

        MatProp(dictionary &transportProperties);
        MatProp(const Time& runTime, const dynamicFvMesh& mesh);

        scalar Sch() const
        {
            return Sch_;
        }

        scalar Pr() const
        {
            return Pr_;
        }

        scalar psat() const
        {
            return psat_;
        }

        scalar cp1() const
        {
            return cp1_;
        }

        scalar cp2() const
        {
            return cp2_;
        }

        scalar k1() const
        {
            return k1_;
        }

        scalar k2() const
        {
            return k2_;
        }

        scalar rho1() const
        {
            return rho1_;
        }

        scalar rho2() const
        {
            return rho2_;
        }

        virtual ~MatProp(){;}
};

}//Foam
#endif


