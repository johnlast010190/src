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
    (c) 2017 Esi Ltd.
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "turbulentTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBaseTurbulenceModel
(
    geometricOneField,
    geometricOneField,
    incompressibleTurbulenceModel,
    IncompressibleTurbulenceModel,
    transportModel
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "laminar/Stokes/Stokes.H"
makeLaminarModel(Stokes);

#include "laminar/Maxwell/Maxwell.H"
makeLaminarModel(Maxwell);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "RAS/SpalartAllmaras/SpalartAllmaras.H"
makeRASModel(SpalartAllmaras);

#include "RAS/EdwardsSA/EdwardsSA.H"
makeRASModel(EdwardsSA);

#include "RAS/kEpsilon/kEpsilon.H"
makeRASModel(kEpsilon);

#include "RAS/RNGkEpsilon/RNGkEpsilon.H"
makeRASModel(RNGkEpsilon);

#include "RAS/realizableKE/realizableKE.H"
makeRASModel(realizableKE);

#include "RAS/LaunderSharmaKE/LaunderSharmaKE.H"
makeRASModel(LaunderSharmaKE);

#include "RAS/kOmega/kOmega.H"
makeRASModel(kOmega);

#include "RAS/kOmegaSST/kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "RAS/buoyantKOmegaSST/buoyantKOmegaSST.H"
makeRASModel(buoyantKOmegaSST);

#include "RAS/kOmegaSSTSAS/kOmegaSSTSAS.H"
makeRASModel(kOmegaSSTSAS);

#include "RAS/kOmegaSSTLM/kOmegaSSTLM.H"
makeRASModel(kOmegaSSTLM);

#include "RAS/v2f/v2f.H"
makeRASModel(v2f);

#include "RAS/LRR/LRR.H"
makeRASModel(LRR);

#include "RAS/SSG/SSG.H"
makeRASModel(SSG);

#include "RAS/kkLOmega/kkLOmega.H"
makeRASModel(kkLOmega);

#include "RAS/ellipticBlendingLagKE/ellipticBlendingLagKE.H"
makeRASModel(ellipticBlendingLagKE);

#include "RAS/ellipticBlendingKE/ellipticBlendingKE.H"
makeRASModel(ellipticBlendingKE);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "LES/Smagorinsky/Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "LES/WALE/WALE.H"
makeLESModel(WALE);

#include "LES/kEqn/kEqn.H"
makeLESModel(kEqn);

#include "LES/dynamicKEqn/dynamicKEqn.H"
makeLESModel(dynamicKEqn);

#include "LES/dynamicLagrangian/dynamicLagrangian.H"
makeLESModel(dynamicLagrangian);

#include "DES/SpalartAllmarasDES/SpalartAllmarasDES.H"
makeLESModel(SpalartAllmarasDES);

#include "DES/SpalartAllmarasDDES/SpalartAllmarasDDES.H"
makeLESModel(SpalartAllmarasDDES);

#include "DES/SpalartAllmarasIDDES/SpalartAllmarasIDDES.H"
makeLESModel(SpalartAllmarasIDDES);

#include "LES/DeardorffDiffStress/DeardorffDiffStress.H"
makeLESModel(DeardorffDiffStress);

#include "DES/kOmegaSSTDES/kOmegaSSTDES.H"
makeLESModel(kOmegaSSTDES);

#include "DES/kOmegaSSTDDES/kOmegaSSTDDES.H"
makeLESModel(kOmegaSSTDDES);

#include "DES/kOmegaSSTIDDES/kOmegaSSTIDDES.H"
makeLESModel(kOmegaSSTIDDES);

#include "DES/kOmegaSSTSASDES/kOmegaSSTSASDES.H"
makeLESModel(kOmegaSSTSASDES);

// ************************************************************************* //
