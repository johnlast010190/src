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
    (c) 2019-2020 Esi Ltd.

Description :
  This file declares the objects and functions for the FoamVOF model
  implementation. FoamVOF is a general-purpose multi-fluid module developed
  under the framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/
#ifndef __FOAMVOF_H
#define __FOAMVOF_H

namespace Foam
{
  class BaroCav;
  class VOFModels;

  class FoamVOF
  {
    public:
      enum {NONE,PISO,PIMPLE,SIMPLE};
    protected:
      const dynamicFvMesh& mesh_;
      const Time& runTime_;
      const MatProp &prop_;
      label mfluid_;
      label solveT_;
      label baroCavitate_;//barotropic cavitation model
      label algorithm_;
      label cavmodel_;
      word name_;
      label use_rgh_;
      scalar g0_;
      bool mules_;
    public:
      bool moveMeshOuterCorrectors;
      FoamVOF(const Time& time, const dynamicFvMesh& mesh,const MatProp &prop);

      List<VOFModels*> models;
      scalar cumulativeContErr;
      label pRefCell;
      scalar pRefValue;
      scalar ghRefv;
      bool LTS;

      bool adjustTimeStep;
      scalar maxCo;
      autoPtr<Function1<scalar> > maxCoDataPtr;
      scalar maxDeltaT;
      scalar minDeltaT;
      scalar CoNum;
      scalar meanCoNum;
      scalar maxAlphaCo;
      scalar alphaCoNum;
      scalar meanAlphaCoNum;

      label nAlphaCorr;
      label nAlphaSubCycles;

      bool MULESCorr;
      bool alphaApplyPrevCorr;
      scalar icAlpha;
      bool alphaRestart ;

      autoPtr<immiscibleIncompressibleTwoPhaseMixture> mixture;
      autoPtr<incompressibleTwoPhaseMixture> mixture2;

      tmp<surfaceScalarField> alphaPhi_;
      tmp<volScalarField> trDeltaT;
      // MULES Correction
      tmp<surfaceScalarField> talphaPhiCorr0;

     virtual surfaceScalarField &alphaPhi();
     virtual volScalarField& alpha1();

     virtual void CourantNo();
     virtual void createAlphaFluxes();
     virtual void alphaCourantNoIF();
     virtual void setRDeltaT(pimpleControl &pimple);
     virtual void createTimeControls();

     void createModels();

     virtual volScalarField& alpha2();

     scalar dt()
     {
        return  runTime_.deltaTValue();
     }

     const MatProp &prop() const
     {
        return prop_;
     }

     bool useRhogh()
     {
         return use_rgh_;
     }

     virtual void createTfields(pimpleControl &pimple);
     virtual void createFields(pimpleControl &pimple);

     virtual void alphaControls();
     virtual void readTimeControls();
     virtual void setDeltaT();
     virtual void setInitialDeltaT();

     autoPtr<Function1<vector>> frameAcceleration;

     void initContinuityErrs()
     {
        cumulativeContErr=0;
     }

     label algorithm()
     {
         return algorithm_;
     }

     label mfluid()
     {
        return mfluid_;//number of fluids
     }

    VOFModels *model(label i)
    {
        return models[i];
    }

     bool hasModel();

     bool mulesScheme()
     {
         return mules_;
     }

     label solveT()
     {
        return solveT_;
     }

     label baroCavitate()
     {
         return baroCavitate_;
     }

     const dynamicFvMesh& mesh()
     {
        return  mesh_;
     }

     const Time& runTime()
     {
        return runTime_;
     }

    label cavitation()
    {
        return cavmodel_;
    }

    volScalarField& p()
    {
        return mesh_.lookupObjectRef<volScalarField>("p");
    }

     volScalarField& p_rgh()
     {
        return mesh_.lookupObjectRef<volScalarField>("p_rgh");
     }

     uniformDimensionedScalarField &hRef()
     {
        return mesh_.lookupObjectRef<uniformDimensionedScalarField>("hRef");
     }

     volScalarField &gh()
     {
        return mesh_.lookupObjectRef<volScalarField>("gh");
     }

     volVectorField& U()
     {
        return mesh_.lookupObjectRef<volVectorField>("U");
     }

     surfaceVectorField& Uf()
     {
        return mesh_.lookupObjectRef<surfaceVectorField>("Uf");
     }

    volScalarField& volTrRate()
    {
        return mesh_.lookupObjectRef<volScalarField>("volTrRate");
    }

    volScalarField& vdotP()
    {
        return mesh_.lookupObjectRef<volScalarField>("vdotP");
    }


    surfaceScalarField &alphaPhiUn()
    {
        return mesh_.lookupObjectRef<surfaceScalarField>("alphaPhiUn");
    }

     volScalarField& rho()
     {
        return mesh_.lookupObjectRef<volScalarField>("rho");
     }

     //heat capacity
     volScalarField& cp()
     {
        return mesh_.lookupObjectRef<volScalarField>("cp");
     }

     //thermal conductivity
     volScalarField& K()
     {
        return mesh_.lookupObjectRef<volScalarField>("K");
     }


     volScalarField& T()
     {
        return mesh_.lookupObjectRef<volScalarField>("T");
     }

     surfaceScalarField &phi()
     {
        return mesh_.lookupObjectRef<surfaceScalarField>("phi");
     }


     surfaceScalarField &rhoPhi()
     {
        return mesh_.lookupObjectRef<surfaceScalarField>("rhoPhi");
     }

     surfaceScalarField &ghf();


     uniformDimensionedVectorField g()
     {
        uniformDimensionedVectorField g1(mesh_.lookupObjectRef<uniformDimensionedVectorField> ("g"));
        if (!use_rgh_)
        {
            g1.value()=0.0;
        }
        return g1;
     }

     virtual const dimensionedScalar& rho1();
     virtual const dimensionedScalar& rho2();

     virtual void alphaEqnSubCycle();

     virtual void UEqnSolve
     (
        fvVectorMatrix& UEqn,
        incompressible::turbulenceModel& turbulence,
        fv::options &fvOptions,
        immiscibleIncompressibleTwoPhaseMixture &mixture,
        pimpleControl &pimple
     );

     virtual void pEqn
     (
        fvVectorMatrix& UEqn,
        fv::options &fvOptions,
        immiscibleIncompressibleTwoPhaseMixture &mixture,
        pimpleControl &pimple
     );

     virtual void updateRho();
     virtual void updateCp();
     virtual void updateK();
     virtual void updateRhoPhi();

     virtual void TEqnSolve
     (
        incompressible::turbulenceModel& turbulence
     );

     virtual void correctAlphav();

     virtual void correctPhi(pimpleControl &pimple);
     virtual void continuityErrs();
     void updateAcceleration(autoPtr<Function1<vector>> frameAcceleration);
     virtual void alphaEqn();

     void dynamicMeshUpdate(pimpleControl &pimple);
     virtual void updateFields
     (
         pimpleControl &pimple,
         incompressible::turbulenceModel &turbulence,
         fv::options& fvOptions
     );

    void alphaEqn(volScalarField &alpha);
    bool hasVolSource();
    bool hasMomentumSource();

    bool hasHeatSource();
    bool hasSpecSource();
    void updateSpec();
    void specSource(volScalarField &Su,volScalarField &Sp);
    void volSource(volScalarField &Su,volScalarField &Sp);

    virtual ~FoamVOF();

   };
}//Foam
#endif
