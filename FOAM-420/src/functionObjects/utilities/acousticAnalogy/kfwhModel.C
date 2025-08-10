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
    (c) 2020-2022 Esi Ltd.

Class
    Foam::functionObjects::kfwhModel

Group
    grpUtilityFunctionObjects

Description
    A founction object which implements an extended aero-acoustics model
    as defined in : Journal of Sound and Vibration (1997), 202 (4), 491-509
    A New Boundary Integral Formulation for the Prediction of Sound Radiation
    By P. di Francescantonio

\*---------------------------------------------------------------------------*/

#include "acousticAnalogy/kfwhModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(kfwhModel, 0);
        addToRunTimeSelectionTable
        (
            functionObject,
            kfwhModel,
            dictionary
        );
    }
}


void Foam::functionObjects::kfwhModel::getVolTerms
(
    bool hasOldField
)
{
    if (modelType_ != "full")
    {
        return;
    }

    const tmp<volSymmTensorField> tR(this->R());
    const volSymmTensorField& tau0 = tR.ref();
    Array2d<scalar> tau;
    //tau.setSize(mesh_.C().size(),6);

    if (viscousEffect_)
    {
        tau.setSize(mesh_.C().size(),6);

        forAll(mesh_.C(), ic)
        {
            const symmTensor& ten = tau0[ic];
            scalar tenf[6] = {ten.xx(),ten.xy(),ten.xz(),ten.yy(),ten.yz(),ten.zz()};
            for (label i=0;i<3;i++)
            {
                for (label j=0;j<3;j++)
                {
                    const label ij0 = 3*i+j;
                    const label i0 = tensorIdx[ij0];
                    tau(ic,i0) = tenf[i0];
                }
            }
        }
    }

    tmp<volScalarField> ptmp = p();
    const volScalarField& pp = ptmp();

    tmp<volVectorField> tU(this->U());
    const volVectorField& UU = tU();
    const tmp<volScalarField> tke = turbEng();

    const scalarField& volumes = mesh_.V();
    scalar dtin = 1.0/dt_;

    label ncells = outCells_.size();

    if (ncells == 0)
    {
        return;
    }

    Tij_.setSize(ncells, 6);

    Tij_ = 0.0;

    scalarField dens = rho();
    for (label ic=0; ic<ncells;ic++)
    {
        if (ic >= pp.size())
        {
            break;
        }
        label ival = outCells_[ic];
        if (ival == 0)
        {
            continue;
        }

        scalar ppp = pp[ic];
        scalar rhoCell = dens[ic];
        const vector& vc = UU[ic];

        for (label i=0;i<3;i++)
        {
            for (label j=0;j<3;j++)
            {
                scalar tauij = 0.0;
                if (viscousEffect_)
                {
                    tauij = tau(ic,ij(i,j));
                }

                scalar pij =
                    (ppp - pref_)*delta(i,j)*presCoef_
                  + rhoCell*vc[i]*vc[j]*veloCoef_;

                scalar tij =
                    pij - 2.0/3.0*rhoCell*tke()[ic]*tkeCoef_*delta(i,j) - tauij;
                Tij_(ic,ij(i,j)) = tij;
            }
        }
    }

    const vectorField& cf = mesh_.C();
    if (!hasOfield_)
    {
        ovterm1_.resize(ncells);
        ovterm2_.resize(ncells);
        dvfield10_.resize(ncells);

        for (label ic=0;ic<ncells;ic++)
        {
            ovterm1_[ic].resize(observerSets_.size());
            ovterm2_[ic].resize(observerSets_.size());
            dvfield10_[ic].resize(observerSets_.size());
            for (label iset=0;iset<observerSets_.size();iset++)
            {
                ovterm1_[ic][iset].resize(observerSets_[iset].size());
                ovterm2_[ic][iset].resize(observerSets_[iset].size());
                dvfield10_[ic][iset].resize(observerSets_[iset].size());
            }
        }
    }

    if (!hasOldField)
    {
        ovterm1_ = 0.0;
        ovterm2_ = 0.0;
    }

    scalar coeff = 1.0/(4.0 * Foam::constant::mathematical::pi);
    scalar currentTime = time_.time().value();

    forAll(observerSets_, obsSetI)
    {
        observer& obsSet = observerSets_[obsSetI];
        const vectorField& pos = obsSet.positions();
        List<scalarList>&  pNonlinearSet = pNonlinear_[obsSetI];
        scalarField& pNonlinearObsCur = obsSet.pNonlinear();

        forAll(pos, obsI)
        {
            scalarList& pnlObs = pNonlinearSet[obsI];
            for (label ic=0; ic<ncells;ic++)
            {
                label ival=outCells_[ic];
                if (ival==0)
                {
                    continue;
                }

                point rvect = pos[obsI] - cf[ic];
                scalar r = mag(rvect);
                if (r<powf(volumes[ic],0.333)*2)
                {
                    continue;
                }

                scalar timeDelay = r/cRef_;
                scalar advancedTime = currentTime + timeDelay;

                scalar trr = Trr(ic,rvect);
                scalar tii = Tii(ic);
                scalar term1 = trr*volumes[ic]/r;
                scalar term2 = (3*trr-tii)*volumes[ic]/(r*r);

                scalar term3 = term3Coef_*(3*trr-tii)*volumes[ic]/(r*r*r);

                scalar vfield1 = term1/(cRef_*cRef_);
                scalar vfield2 = term2/cRef_;
                scalar vfield3 = term3;

                scalar dvfield1 = 0.0;
                scalar dvfield2 = 0.0;
                scalar ddvfield1 = 0.0;
                if (hasOldField)
                {
                    dvfield1 = dtin*(vfield1 - ovterm1_(ic, obsSetI, obsI));
                    dvfield2 = dtin*(vfield2 - ovterm2_(ic, obsSetI, obsI));
                }
                if (hasOOfield_)
                {
                    ddvfield1 = dtin*(dvfield1 - dvfield10_(ic, obsSetI, obsI));
                }
                scalar pnlTerm=ddvfield1+dvfield2+vfield3;
                if (advancedTime < endSamplingTime_)
                {
                    label previousTimeStep =
                            floor
                            (
                                (advancedTime-startSamplingTime_)
                              / dtAco_
                            );

                    scalar previousTime = startSamplingTime_
                            + previousTimeStep*dtAco_;

                    scalar weight =
                        1.0 - (advancedTime - previousTime)/dtAco_;
                    pnlObs[previousTimeStep] += weight*pnlTerm;
                    if (previousTimeStep + 1 < pnlObs.size())
                    {
                         pnlObs[previousTimeStep+1] += (1-weight)*pnlTerm;
                    }
                }

                ovterm1_(ic, obsSetI, obsI) = vfield1;
                ovterm2_(ic, obsSetI, obsI) = vfield2;
                dvfield10_(ic,obsSetI,obsI) = dvfield1;
            }


            label currentTimeStep =
                floor
                (
                    (time_.time().value() - startSamplingTime_)/dtAco_
                );

            scalar pnl = pnlObs[currentTimeStep];
            reduce(pnl, sumOp<scalar>());
            pNonlinearObsCur[obsI] = coeff*pnl;
        }
    }

    return;
}


Foam::functionObjects::kfwhModel::kfwhModel
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    FfowcsWilliamsHawkings
    (
        name,
        runTime,
        dict,
        "KFWH"
    )
{
    veloCoef_ = 1.0; //whether to consider velocity in the volume integral
    presCoef_ = 1.0;
    tkeCoef_ = 0.0;
    term3Coef_ = 1.0;

    phaseNo_ = 1;
    boundSurfaceFile_ = "boundSurface.stl";
    pref_ = 0.0;
    //modelType_="simplified"; //or simplified

    meshChange_ = false;
    prefOpt_ = "mean";
    icPref_ = 0;

    hpcMode_ = 1;
    read(dict);

    pThickness_.setSize(observerSets_.size());
    pLoad_.setSize(observerSets_.size());
    pNonlinear_.setSize(observerSets_.size());

    getSampleBox();
    getFaceNormals();

    Info<< "build surface KNN..." << endl;
    buildSurfaceKnn();

    getCellList();

    KNN cellMapKnn;
    buildCellMapKnn(cellMapKnn);

    getFaceCellMap(cellMapKnn);

    insolid_.resize(faceToCell_.size(), 0);

    L_.setSize(faceToCell_.size(), 6);
    L0_.setSize(faceToCell_.size(), 6);
    L00_.setSize(faceToCell_.size(), 6);
    dotL_.setSize(faceToCell_.size(), 6);

    Uf_.resize(faceToCell_.size());
    Uf0_.resize(faceToCell_.size());
    Uf00_.resize(faceToCell_.size());

    if (modelType_ == "full")
    {
        outCells_.resize(mesh_.C().size(), 0);

        forAll(mesh_.C(),i)
        {
            const point& pt = mesh_.C()[i];
            if (!inSurfObject(pt))
            {
                outCells_[i] = 1;
            }
        }
    }

    dt_ = mesh_.time().deltaTValue();
    dto_ = dt_;

    hasOfield_ = false;
    hasOOfield_ = false;

    if (prefOpt_ == "fixloc")
    {
        std::vector<label> nb;
        std::vector<scalar> dist;
        std::vector<scalar> qvect(3);
        qvect[0] = prefLoc_.x();
        qvect[1] = prefLoc_.y();
        qvect[2] = prefLoc_.z();

        cellMapKnn.search
        (
            qvect,
            nb,
            dist
        );
        icPref_ = nb[0];
    }
    setSurfaceState();
}



void Foam::functionObjects::kfwhModel::setSurfaceState()
{
    insolid_ = 0;
    if (sobject_.type_ == "unknown")
    {
        return;
    }
    const vectorField& cf = Cf();
    forAll(cf, faceI)
    {
        const point& pt = cf[faceI];
        if (sobject_.inObject(pt))
        {
            insolid_[faceI] = 1;
        }
    }

    label isum = 0;
    for (label i=0;i<insolid_.size();i++)
    {
        isum += insolid_[i];
    }
    Info<< "Number of faces in solid: " << isum
        << ",  total faces: " << insolid_.size() << endl;

}


void Foam::functionObjects::kfwhModel::initialise()
{
    FfowcsWilliamsHawkings::initialise();

    scalar totalNumberOfTimeSteps = endSamplingTime_-startSamplingTime_;

    totalNumberOfTimeSteps =
        min
        (
            floor(totalNumberOfTimeSteps/dtAco_) + 1,
            scalar(1E7)
        );

    label sizeTime = floor(totalNumberOfTimeSteps);

    label sizeThk = sizeTime;

    forAll(observerSets_, obsI)
    {
        pThickness_[obsI].setSize
        (
            observerSets_[obsI].size(),
            scalarField(sizeThk, 0.0)
        );
        pLoad_[obsI].setSize
        (
            observerSets_[obsI].size(),
            scalarField(sizeThk, 0.)
        );
        pNonlinear_[obsI].setSize
        (
            observerSets_[obsI].size(),
            scalarField(sizeThk, 0.)
        );
    }
}


void Foam::functionObjects::kfwhModel::writeFfowcsWilliamsHawkings()
{
    Log << "Sum of instantaneous acoustic pressure: " << endl;
    forAll(observerSets_, obsSetI)
    {
        observer& obs = observerSets_[obsSetI];
        scalar pPrimeSum = gSum(obs.pPrime());
        Log << "   -" << obs.name() << ": " << pPrimeSum << " [Pa]" << endl;
    }

    // File output
    //file() << obr_.time().timeName() << tab << setw(1) << "   ";
    label currentTimeStep = floor
    (
        (time_.time().value() - startSamplingTime_) / dtAco_
    );
    scalar stime=currentTimeStep*dtAco_;

   // file() << obr_.time().timeName() << tab << setw(1) << "   ";
    file() <<stime<< tab << setw(1) << "   ";

    forAll(observerSets_, obsSetI)
    {
        observer& obs = observerSets_[obsSetI];
        if (obs.write())
        {
            file() << obs.name() << "    ";

            if (writeComponent_)
            {
                const scalarField& pPrime = obs.pPrime();
                const scalarField& pthick = obs.pThickness();
                const scalarField& pload = obs.pLoad();
                const scalarField& pnonlinear = obs.pNonlinear();
                forAll(obs.positions(), obsI)
                {
                    file() << pPrime[obsI] ;
                    if (writeComponent_)
                    {
                        file()<< "( "<< pthick[obsI]
                        <<"  "<<pload[obsI]
                        <<"  "<<pnonlinear[obsI]<<")"
                        << "   ";
                    }
                }
            }
            else
            {
                const scalarField& pPrime = obs.pPrime();
                forAll(obs.positions(), obsI)
                {
                    file() << pPrime[obsI] << "   ";
                }
            }
        }
    }
    file() << endl;
}


void Foam::functionObjects::kfwhModel::getCellList()
{
    Info<< "Getting cell list..." << endl;
    vectorField searchZoneXyz;
    labelList cellIds;
    labelList procIds;
    label myid = Pstream::myProcNo();
    forAll(mesh_.C(),ic)
    {
        const point& pt = mesh_.C()[ic];
        if (inSearchZone(pt))
        {
            searchZoneXyz.append(pt);
            cellIds.append(ic);
            procIds.append(myid);
        }
    }

    if (!Pstream::parRun())
    {
        searchZoneXyz_ = searchZoneXyz;
        cellIds_ = cellIds;
        procIds_ = procIds;
    }
    else
    {
        //coordinates
        List<vectorField> posLoc(Pstream::nProcs());
        posLoc[myid] = searchZoneXyz;
        Pstream::allGatherList(posLoc); //make it available for other processors

        searchZoneXyz_ = ListListOps::combine<vectorField>(posLoc, accessOp<vectorField>());
        //cellId
        List<labelList> cellIdLoc(Pstream::nProcs());
        cellIdLoc[myid] = cellIds;
        Pstream::allGatherList(cellIdLoc); //make it available for other processors

        cellIds_ = ListListOps::combine<labelList>(cellIdLoc, accessOp<labelList>());
        //processor id
        List<labelList> procIdLoc(Pstream::nProcs());
        procIdLoc[myid] = procIds;
        Pstream::allGatherList(procIdLoc); //make it available for other processors

        procIds_ = ListListOps::combine<labelList>(procIdLoc, accessOp<labelList>());
    }
    Info<< "cells in search zone:" << searchZoneXyz_.size() << " "
        << cellIds_.size() << " " << procIds_.size() << endl;
}


void Foam::functionObjects::kfwhModel::getFaceCellMap(KNN &knn)
{
    const vectorField& cf = Cf();
    faceToCell_.resize(cf.size());

    forAll(cf,j)
    {
        const point& pt = cf[j];
        std::vector<label> nb;
        std::vector<scalar> dist;
        std::vector<scalar> qvect(3);
        qvect[0] = pt.x();
        qvect[1] = pt.y();
        qvect[2] = pt.z();

        knn.search
        (
            qvect,
            nb,
            dist
        );
        label i0 = nb[0];
        faceToCell_[j] = i0;
    }
    Info<< "face to cell map obtained." << endl;
}


void Foam::functionObjects::kfwhModel::getFaceNormals()
{
    Info<< "getting surface normals...." << endl;

    const vectorField &cf=Cf();
    point vmean(0,0,0);
    forAll(cf,j)
    {
        vmean += cf[j];
    }
    vmean /= float(Cf().size());

    faceNormals_ = Sf()/magSf();
    forAll(cf, j)
    {
        //point nsf=Sf()[j]/magSf()[j];
        point nsf = faceNormals_[j];
        point dsf = cf[j] - vmean;
        scalar d = nsf&dsf;
        if (d<0)
        {
            nsf *= (-1);
            faceNormals_[j] = nsf;
        }
    }
}


void Foam::functionObjects::kfwhModel::buildSurfaceKnn()
{
    Info<< "start build surface knn..." << endl;
    const vectorField& cf = Cf();
    std::vector<std::vector<scalar>> trainData;

    forAll(cf, j)
    {
        const point pt = cf[j];
        std::vector<scalar> pts(3);
        pts[0] = pt.x();
        pts[1] = pt.y();
        pts[2] = pt.z();
        trainData.push_back(pts);
    }
    surfaceKnn_.build(trainData);
    Info<< "knn build." << endl;
}


void Foam::functionObjects::kfwhModel::buildCellMapKnn(KNN &knn)
{
    Info<< "start build cellMap knn..." << endl;

    std::vector<std::vector<scalar>> trainData;
    if (!Pstream::parRun()||hpcMode_==1)
    {
        forAll(searchZoneXyz_,j)
        {
            const point& pt = searchZoneXyz_[j];
            std::vector<scalar> pts(3);
            pts[0] = pt.x();
            pts[1] = pt.y();
            pts[2] = pt.z();
            trainData.push_back(pts);
        }
    }
    else
    {
        List<vectorField> posLoc(Pstream::nProcs());
        label myid = Pstream::myProcNo();
        posLoc[myid] = mesh_.C();
        Pstream::allGatherList(posLoc); //make it available for other processors

        vectorField posGlob = ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());

        forAll(posGlob, i)
        {
            const point& pt = posGlob[i];
            std::vector<scalar> pts(3);
            pts[0] = pt.x();
            pts[1] = pt.y();
            pts[2] = pt.z();
            trainData.push_back(pts);
        }
    }

    knn.build(trainData);

    Info<< "cellMap knn build." << endl;
}


Foam::word Foam::functionObjects::kfwhModel::surfaceFile()
{
    word root = cwd();
    word dname = root + "/constant/triSurface/" + boundSurfaceFile_;
    return dname;
}


void Foam::functionObjects::kfwhModel::setTensorField()
{
    Array2d<scalar> tau;
    tau.setSize(mesh_.C().size(), 6);
    tau = 0.0;

    const tmp<volSymmTensorField> tR(this->R());
    const volSymmTensorField& tau0 = tR.ref();
    if (viscousEffect_)
    {
        forAll(mesh_.C(),ic)
        {
            const symmTensor& ten = tau0[ic];
            scalar tenf[6] =
                {ten.xx(), ten.xy(), ten.xz(), ten.yy(), ten.yz(), ten.zz()};
            for (label i=0;i<3;i++)
            {
                for (label j=0;j<3;j++)
                {
                    label ij0 = 3*i+j;
                    label i0 = tensorIdx[ij0];
                    tau(ic,i0) = tenf[i0];
                }
            }
        }
    }

    tmp<volScalarField> ptmp = p();
    const volScalarField& pp = ptmp();

   // scalarField pp(p());

    tmp<volVectorField> tU(this->U());
    const volVectorField& UU = tU();

    Info<< "setting tensor field..." << tau.size() << endl;

    L_ = 0.0;

    label myid = Pstream::myProcNo();
    vector vzero(0,0,0);
    Uf_ = vzero;
    Info<< "faceToCell size:" << faceToCell_.size() << endl;

    scalarField dens=rho();
    if (!Pstream::parRun()||hpcMode_==1)
    {
        forAll(faceToCell_,jf)
        {
            label idx = faceToCell_[jf];
            label ic = cellIds_[idx];
            label procId = procIds_[idx];

            bool toAdd = false;
            if (!Pstream::parRun() || procId == myid)
            {
                toAdd = true;
            }

            if (toAdd)
            {
                scalar ppp = pp[ic];
                const vector& vc = UU.primitiveField()[ic];
                scalar rhoCoef = dens[ic]/rhoRef_;

                for (label i=0;i<3;i++)
                {
                    for (label j=0;j<3;j++)
                    {
                        scalar tauij = tau(ic,ij(i,j));
                        scalar pij = (ppp-pref_)*delta(i,j) - tauij;
                        L_(jf,ij(i,j)) = pij+dens[ic]*vc[i]*vc[j];
                    }
                }
                Uf_[jf] = rhoCoef*vc;
            }
        }//forAll
    }
    else
    {
        List<scalarField> locP(Pstream::nProcs());
        locP[myid] = pp;
        Pstream::allGatherList(locP);
        scalarField pGlobal=ListListOps::combine<scalarField>(locP,accessOp<scalarField>());
        List<vectorField> locU(Pstream::nProcs());
        locU[myid] = UU;
        Pstream::allGatherList(locU);
        scalar coefs = 0;
        if (myid == 0)
        {
            coefs=1;
        }

        vectorField UGlobal = ListListOps::combine<vectorField>(locU,accessOp<vectorField>());
        forAll(faceToCell_, jf)
        {
            label ic = faceToCell_[jf];
            scalar ppp = pGlobal[ic];
            const vector& vc = UGlobal[ic];
            for (label i=0;i<3;i++)
            {
                for (label j=0;j<3;j++)
                {
                    scalar pij = (ppp - pref_)*delta(i,j);
                    L_(jf,ij(i,j)) = coefs*(pij + rhoRef_*vc[i]*vc[j]);
                }
            }

            Uf_[jf] = coefs*vc;
        }
    }
    Info<< "tensor field set." << endl;
    return;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::kfwhModel::turbEng()
{
    typedef incompressible::turbulenceModel icoTurbModel;
    typedef compressible::turbulenceModel cmpTurbModel;

    if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel &turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);
        return turb.k();
    }
    else if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const compressible::turbulenceModel &turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);
        return turb.k();
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


const Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::kfwhModel::R()
{
    typedef incompressible::turbulenceModel icoTurbModel;
    typedef compressible::turbulenceModel cmpTurbModel;

    if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel &turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return turb.devReff();
    }
    else if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel &turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        viscousEffect_ = false;
        return volSymmTensorField::null();
    }
}


void Foam::functionObjects::kfwhModel::getSampleBox()
{
    const vectorField& cf = Cf();
    boundBox bbox(cf);

    point dpt = bbox.max() - bbox.min();
    const scalar eps = 0.05;
    const scalar xmin = bbox.min().x() - dpt.x()*eps;
    const scalar xmax = bbox.max().x() + dpt.x()*eps;
    const scalar ymin = bbox.min().y() - dpt.y()*eps;
    const scalar ymax = bbox.max().y() + dpt.y()*eps;
    const scalar zmin = bbox.min().z() - dpt.z()*eps;
    const scalar zmax = bbox.max().z() + dpt.z()*eps;
    point pmin(xmin, ymin, zmin);
    point pmax(xmax, ymax, zmax);
    samplebox_ = new boundBox(pmin,pmax);
    Info<< "bound box: " << samplebox_().min() << "," << samplebox_().max() << endl;
    return;
}


bool Foam::functionObjects::kfwhModel::read(const dictionary& dict)
{
    term3Coef_ = 1.0;
    presCoef_ = 1.0;
    veloCoef_ = 1.0;
    tkeCoef_ = 1.0;

    boundSurfaceFile_=dict.lookupOrDefault<word>("boundSurfaceFile", "boundSurface.stl");

    meshChange_ = dict.lookupOrDefault<bool>("meshChange", false);
    pref_ = dict.lookupOrDefault<scalar>("pref", 0.0);
    viscousEffect_ = dict.lookupOrDefault<bool>("visStress", false);
    phaseNo_ = dict.lookupOrDefault<label>("phaseNo", 1);

    hpcMode_ = dict.lookupOrDefault<label>("hpcMode", 1);

    Info<< "HPC mode..." << hpcMode_ << endl;

    //whether to include pressure  in the quadrpole  source
    bool qpresTerm = dict.lookupOrDefault<bool>("qPresSource", true);

    if (!qpresTerm)
    {
         presCoef_ = 0.0;
    }

    bool velTerm = dict.lookupOrDefault<bool>("qVelSource", true);
    if (!velTerm)
    {
        veloCoef_ = 0.0;
    }

    if (velTerm)
    {
        tkeCoef_ = 0;
    }

    //mean,fixed,fixloc
    prefOpt_ = dict.lookupOrDefault<word>("prefOption", "fixed");

    bool sobject = dict.lookupOrDefault<bool>("constraint", false);

    if (prefOpt_=="fixloc")
    {
        Info<< "reference pressure uses fixed location." << endl;
        prefLoc_ = dict.lookupOrDefault<point>("prefLoc", prefLoc_);
    }

    word bfile = surfaceFile();
    surfaces_ = new triSurface(obr_,bfile);
    // read solid object
    if (sobject)
    {
        const dictionary& objDict = dict.subDict("solidObject");
        word name = objDict.lookupOrDefault<word>("name", "none");
        sobject_.name_ = name;
        word type = objDict.lookupOrDefault<word>("type", "unknown");
        sobject_.type_ = type;

        vectorField pts;
        sobject_.points_ = objDict.lookupOrDefault<vectorField>("points", pts);
        sobject_.r_ = objDict.lookupOrDefault<scalar>("r", 1.0);

        if (sobject_.points_.size() < 1)
        {
            FatalErrorInFunction
                << "Number of points must be at least 1 to define a solidObject "
                << exit(FatalError);
        }

        label gridCells = objDict.lookupOrDefault<label> ("gridCells",40000);

        if (type == "smBoundBox")
        {
            if (sobject_.points_.size()<=1)
            {
                FatalErrorInFunction
                << "Minimum and maximum vertices need to be defined for smBoundBox"
                << exit(FatalError);
            }

            const point& pmin = sobject_.points_[0];
            const point& pmax = sobject_.points_[1];
            sobject_.solidBox_.reset
            (
                new smBoundBox
                (
                    pmin,
                    pmax,
                    "smBoundBox",
                    1
                )
            );

            sobject_.solidBox_().getCellList(&mesh_,gridCells);
        }
    }

    if (modelType_ == "full")
    {
        writeComponent_ = true;
    }

    if (term3Coef_ < 0.1)
    {
        veloCoef_ = 1;
        presCoef_ = 1;
    }

    return true;
}


Foam::scalarField Foam::functionObjects::kfwhModel::rho()
{
    label sz = mesh_.V().size();
    scalarField rhofld(sz, rhoRef_);
    if (phaseNo_ == 2)
    {
        const volScalarField& dens = lookupObject<volScalarField>("rho");
        forAll(dens.primitiveField(),i)
        {
            rhofld[i] = dens.primitiveField()[i];
        }
    }
    return rhofld;
}


void Foam::functionObjects::kfwhModel::calculate()
{
    Info<< "doing calc..." << endl;
    //scalar tcurrent= time_.time().value();
    //if (tcurrent < startSamplingTime_ || tcurrent > endSamplingTime_)
    //{
    //    return;
    //}

    initialise();
    //reference pressure
    label myid = Pstream::myProcNo();

    tmp<volScalarField> tpres(this->p());
    const volScalarField& pres = tpres();

    dt_ = mesh_.time().deltaTValue();

    scalarField pvol(pres.primitiveField()*mesh_.V());

    scalar pSum = gSum(pvol);

    labelList cellNo(Pstream::nProcs());

    label ncells = mesh_.C().size();
    cellNo[myid] = ncells;

    //label totCells=gSum(cellNo);
    //scalar pmean=pSum/float(totCells);
    scalar pmean = pSum/gSum(mesh_.V());
    if (prefOpt_ == "mean")
    {
        pref_ = pmean;
    }

    if (prefOpt_ == "fixloc")
    {
        pref_ = pres[icPref_];
        Info<< "reference pressure cell and pref: "
            << icPref_ << " " << pref_ << endl;
    }

    if (meshChange_)
    {
        KNN cellMapKnn;
        buildCellMapKnn(cellMapKnn);
        getFaceCellMap(cellMapKnn);
    }

    Info<< "mean pressure:" << ncells << " " << pmean << endl;

    setTensorField();
    //integrations
    scalar coeff = 1.0/(4.0*Foam::constant::mathematical::pi);

    vectorField dotUf(Uf_);

    dotUf = 0.0;
    dotL_ = 0.0;

    //time deritives
    scalar coefft   = 1.0/dt_ + 1.0/(dt_ + dto_);
    scalar coefft00 = dt_/(dto_*(dt_ + dto_));
    scalar coefft0  = coefft + coefft00;

   // scalar dtin=1.0/dt_;
    if (modelType_ == "full")
    {
        Info<< "get volume integrals..." << endl;
        getVolTerms(hasOfield_);
    }

    Info<< "Assemble results...." << endl;
    if (hasOOfield_)
    {
        for (label i=0;i<dotL_.dim1();i++)
        {
            for (label j=0;j<dotL_.dim2();j++)
            {
                dotL_(i,j) =
                    coefft*L_(i,j) - coefft0*L0_(i,j) + coefft00*L00_(i,j);
            }
        }
        dotUf =
            coefft*Uf_ - coefft0*Uf0_ + coefft00*Uf00_;
    }

    const scalar Cinv = 1.0/cRef_;
    scalar currentTime = time_.time().value();


    forAll(observerSets_, obsSetI)
    {
        List<scalarList>& pPrimeSet = pPrime_[obsSetI];
        List<scalarList>& pThicknessSet = pThickness_[obsSetI];
        List<scalarList>& pLoadSet = pLoad_[obsSetI];

        observer& obsSet = observerSets_[obsSetI];
        const vectorField& pos = obsSet.positions();

        scalarField& pPrimeObsCur = obsSet.pPrime();
        scalarField& pThicknessObsCur = obsSet.pThickness();
        scalarField& pLoadObsCur = obsSet.pLoad();
        //scalarField& pNonlinearObsCur = obsSet.pNonlinear();

        forAll(pos, obsI)
        {
            scalar pPrime = 0.0;
            scalar pthickness = 0.0;
            scalar pload = 0.0;

            scalarList& pPrimeObs = pPrimeSet[obsI];
            scalarList& pThicknessObs = pThicknessSet[obsI];
            scalarList& pLoadObs = pLoadSet[obsI];

            forAll(Cf(), jf)
            {
                if (insolid_[jf] == 1)
                {
                    continue;
                }

                const point& cf = Cf()[jf];
                point rvect = pos[obsI] - cf;
                scalar r = mag(rvect);
                scalar timeDelay = r/cRef_;
                scalar advancedTime = currentTime + timeDelay;

                if (advancedTime < endSamplingTime_)
                {
                    scalar term1 = rhoRef_*dotUf[jf]&Sf()[jf]/r;
                    pthickness = term1;

                    scalar term2 = Cinv*magSf()[jf]*dFiri(jf,rvect)/r;
                    scalar term3 = magSf()[jf]*Firi(jf,rvect)/(r*r);
                    pload = term2 + term3;

                    pPrime = term1 + term2 + term3;

                    label previousTimeStep =
                        floor
                        (
                            (advancedTime - startSamplingTime_)/dtAco_
                        );

                    const scalar previousTime =
                        startSamplingTime_ + previousTimeStep*dtAco_;

                    const scalar weight =
                        1.0 - (advancedTime - previousTime)/dtAco_;

                    pPrimeObs[previousTimeStep] += weight*pPrime;
                    pThicknessObs[previousTimeStep] += weight*pthickness;
                    pLoadObs[previousTimeStep] += weight*pload;

                    if (previousTimeStep + 1 < pPrimeObs.size())
                    {
                        pPrimeObs[previousTimeStep+1] += (1 - weight)*pPrime;
                        pThicknessObs[previousTimeStep+1] += (1 - weight)*pthickness;
                        pLoadObs[previousTimeStep+1] += (1 - weight)*pload;
                    }
                }
            }//for all cf

            label currentTimeStep = floor
            (
                (time_.time().value() - startSamplingTime_) / dtAco_
            );

            pPrime = pPrimeObs[currentTimeStep];

            pthickness = pThicknessObs[currentTimeStep];
            pload = pLoadObs[currentTimeStep];
            reduce(pPrime, sumOp<scalar>());
            reduce(pthickness, sumOp<scalar>());
            reduce(pload, sumOp<scalar>());
            pPrimeObsCur[obsI] = coeff*pPrime;

            pThicknessObsCur[obsI] = coeff*pthickness;
            pLoadObsCur[obsI] = coeff*pload;
        }
    }

    // set old values
    if (hasOfield_)
    {
        L00_ = L0_;
        Uf00_ = Uf0_;
        hasOOfield_ = true;
    }

    L0_ = L_;
    Uf0_ = Uf_;

    hasOfield_ = true;
    dto_ = dt_;
}


// ************************************************************************* //
