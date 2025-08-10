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

  This file declares the objects and functions for cavitation models used
  in FoamVOF model implementation. FoamVOF is a general-purpose multi-fluid
  module developed under the framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/

#ifndef __CAVITATE_H
#define __CAVITATE_H
#include "vofmodels.h"

namespace Foam
{
class FoamVOF ;
class Cavitate
:
   public VOFModels
{
    protected:
       scalar nseed_;
       scalar R0_;
       scalar fvmin_;//minimum vapour fraction
       scalar fvmax_;
       word vapName_;
       word liqName_;
       scalar rhol_;
       scalar rhov_;
    public:
       Cavitate
       (
           const MatProp &prop,
           const FoamVOF &vof
       );

       void volTransRate(volScalarField &vtr);
       const volScalarField& vapour();

       scalar R0()
       {
           return R0_;
       }

       scalar fvconst()
       {
         return 4.0*constant::mathematical::pi*nseed_/3.0;
       }

       void updateSpec()
       {
          volTransRate(vof().volTrRate());
       }

       label massTransfer()
       {
          return 1;
       }

       label specTransfer()
       {
          return 1;
       }

       scalar sign(scalar f)
       {
           if (f>=0)
           {
             return 1.0;
           }
           else
           {
             return -1.0;
           }
        }

        void specSource(const word& scname, volScalarField &Su,volScalarField &Sp);

        virtual ~Cavitate(){;}

};


}// End namespace Foam

#endif

