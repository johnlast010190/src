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
  This file declares the variables and functions for the MulesAlgo class in
  FoamVOF model implementation. FoamVOF is a general-purpose multi-fluid
  module developed under the framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/
#ifndef MulesAlgo_H
#define MulesAlgo_H

namespace Foam
{
class MulesAlgo : public FoamVOF
{
    public:
        MulesAlgo
        (
            const Time& time,
            const dynamicFvMesh& mesh,
            const MatProp &prop
        );

        void updateFields
        (
            pimpleControl &pimple,
            incompressible::turbulenceModel &turbulence,
            fv::options& fvOptions
        );

        void createFields
        (
            pimpleControl &pimple
        );

    virtual ~MulesAlgo(){;}
};

//******************************************//
} //Foam
#endif






