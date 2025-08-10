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
  This file declares the base class for physical models in FoamVOF
  implementation. FoamVOF is a general-purpose multi-fluid module developed
  under the framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/
#ifndef __VOFMODELS_H
#define __VOFMODELS_H

namespace Foam
{

class VOFModels
{
    protected:
       const MatProp &prop_;
       const FoamVOF &vof_;
       word modelName_;
    public:
       VOFModels
       (
          const MatProp &prop,
          const FoamVOF &vof,
          const word &name="none"
       );

       const MatProp& prop()
       {
           return prop_;
       }

       const word& name()
       {
           return  modelName_;
       }

       const Time &runTime()
       {
           return const_cast<FoamVOF &>(vof()).runTime();
       }

       const dynamicFvMesh& mesh()
       {
           return const_cast<FoamVOF &>(vof()).mesh();
       }

       const FoamVOF& vof() const
       {
           return vof_;
       }

      FoamVOF& vof()
      {
          return const_cast<FoamVOF&>(vof_);
      }

      virtual label massTransfer()
      {
          return 0;
      }

      virtual label heatTransfer()
      {
          return 0;
      }

      virtual label specTransfer()
      {
          return 0;
      }

      virtual label momentumTransfer()
      {
          return 0;
      }
      virtual void updateSpec(){;}

      virtual void volSource(volScalarField &Su,volScalarField &Sp);
      virtual void heatSource(volScalarField &Su,volScalarField &Sp);
      virtual void specSource(const word& scname, volScalarField &Su,volScalarField &Sp);
      virtual void velSource(volVectorField &Su,volScalarField &Sp);

      virtual ~VOFModels(){;}
};

}// End namespace Foam

#endif

