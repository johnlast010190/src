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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/meshMetric/meshMetric.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::meshMetric::sphericityCalculation(const label& pI)
{
    const labelList& cells = mesh_.pointCells()[pI];

    spher_.setSize(cells.size(), 0);

    spher_ *= 0;
    nonSphericity_ = spher_;
    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        spher_[cI] = sphericityValue(cL);
    }
}

Foam::scalar Foam::meshMetric::sphericityValue(const label& cL)
{
    const scalar& pi = constant::mathematical::pi;
    scalar spher;

    scalar vol = cellVols()[cL];
    scalar surf = cellSurface()[cL];
    scalar numerator = vol*vol;
    scalar denominator = surf*surf*surf;
    spher = 36*pi*numerator/denominator*cellTypeCoef()[cL];

    //add something small to sphericity to limit the derivatives
    //of 1/spher in magnitude
    spher += 0.001;

    return spher;
}

Foam::vector Foam::meshMetric::sphericityDer
(
    const scalar& Vol, const scalar& Surf,
    const vector& d_Vol, const vector& d_Surf
)
{
    vector calc = vector::zero;
    scalar vol = Vol;// + VSMALL;
    const scalar& pi = constant::mathematical::pi;
    scalar coef1 = 2.;
    scalar coef2 = 3.*(vol)/(Surf);
    scalar coef3 = 36*pi*vol/(Surf*Surf*Surf);
    calc = coef3*(coef1*d_Vol - coef2*d_Surf);
    return calc;
}

Foam::tensor Foam::meshMetric::sphericity2ndDer
(
    const scalar& vol, const scalar& surf,
    const vector& dVol, const vector& dSurf,
    const tensor& d_2Sur
)
{
    tensor t = tensor::zero;

    t.xx() = 2*dVol.x()*dVol.x() + 12*dSurf.x()*dSurf.x()/surf/surf*vol*vol
            -6*vol/surf*(dSurf.x()*dVol.x()+dSurf.x()*dVol.x())
            -3*vol*vol/surf*(d_2Sur.xx());

    t.yy() = 2*dVol.y()*dVol.y() + 12*dSurf.y()*dSurf.y()/surf/surf*vol*vol
            -6*vol/surf*(dSurf.y()*dVol.y()+dSurf.y()*dVol.y())
            -3*vol*vol/surf*(d_2Sur.yy());

    t.zz() = 2*dVol.z()*dVol.z() + 12*dSurf.z()*dSurf.z()/surf/surf*vol*vol
            -6*vol/surf*(dSurf.z()*dVol.z()+dSurf.z()*dVol.z())
            -3*vol*vol/surf*(d_2Sur.zz());

    t.xy() = 2*dVol.x()*dVol.y() + 12*dSurf.x()*dSurf.y()/surf/surf*vol*vol
            -6*vol/surf*(dSurf.x()*dVol.y()+dSurf.y()*dVol.x())
            -3*vol*vol/surf*(d_2Sur.xy());

    t.xz() = 2*dVol.x()*dVol.z() + 12*dSurf.x()*dSurf.z()/surf/surf*vol*vol
            -6*vol/surf*(dSurf.x()*dVol.z()+dSurf.z()*dVol.x())
            -3*vol*vol/surf*(d_2Sur.xz());

    t.yz() = 2*dVol.y()*dVol.z() + 12*dSurf.y()*dSurf.z()/surf/surf*vol*vol
            -6*vol/surf*(dSurf.y()*dVol.z()+dSurf.y()*dVol.x())
            -3*vol*vol/surf*(d_2Sur.yz());

    t.yx() = t.xy();
    t.zy() = t.yz();
    t.zx() = t.xz();

    const scalar& pi = constant::mathematical::pi;
    scalar coef = 36*pi/(surf*surf*surf);
    t *= coef;

    return t;
}

void Foam::meshMetric::nonSphericityHessian
(
    tensor& t,
    const vector& der,
    const label& cI
)
{
    label pow = pow_ + 4;
    scalar coef1 = nonSphericity_[cI]/spher_[cI];
    scalar coef2 = coef1/spher_[cI];
    t.xx() = pow*coef2*der.x()*der.x() - coef1*t.xx();
    t.yy() = pow*coef2*der.y()*der.y() - coef1*t.yy();
    t.zz() = pow*coef2*der.z()*der.z() - coef1*t.zz();
    t.xy() = pow*coef2*der.x()*der.y() - coef1*t.xy();
    t.xz() = pow*coef2*der.x()*der.z() - coef1*t.xz();
    t.yz() = pow*coef2*der.y()*der.z() - coef1*t.yz();
    t.yx() = t.xy();
    t.zy() = t.yz();
    t.zx() = t.xz();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::meshMetric::meshMetric
(
    const fvMesh& mesh_,
    pointField& state,
    const dictionary& dict
)
:
    meshGeometry(mesh_, state, dict),
    mesh_(mesh_),
    spher_(0),
    nonSphericity_(0),
    objective_(0),
    derivative_(vector::zero),
    hessian_(tensor::zero),
    minSphericityCell_(-1),
    pow_(2)
{
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //
Foam::meshMetric::~meshMetric()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::scalar Foam::meshMetric::pointValue(const label& pI)
{
    sphericityCalculation(pI);
    objective_ = 0;

    const labelList& cells = mesh_.pointCells()[pI];
    forAll(cells, cI)
    {
        scalar minValue = 1.0e-50;
        if (spher_[cI]<minValue)
        {
            nonSphericity_[cI] = VGREAT;
        }
        else
        {
            nonSphericity_[cI] = spher_[cI];
            for (int i=0; i<pow_-1; i++)
            {
                nonSphericity_[cI] *= spher_[cI];
            }
            nonSphericity_[cI] = 1./nonSphericity_[cI];
        }
    }

    forAll(cells, cI)
    {
        objective_ += nonSphericity_[cI];
    }

    objective_ /= pointCellsSize()[pI];

    return objective_;
}

Foam::tmp<Foam::scalarField> Foam::meshMetric::pointValueGlobal(const labelList& pList)
{
    tmp<scalarField> tobjective(new scalarField(pList.size()));
    scalarField& objective = tobjective.ref();

    for (int I=0; I<pList.size(); I++)
    {
        const label& pI = pList[I];
        objective[I] = pointValue(pI);
    }

    syncTools::syncPointList
    (
        mesh_,
        pList,
        objective,
        plusEqOp<scalar>(),
        scalar(0)
    );

    return tobjective;
}

Foam::vector Foam::meshMetric::displacement(const label& pI)
{
    calculateGradients(pI);

    vector disp = modifiedNewtonsMethod::newtonsStep(derivative_, hessian_);

    limitDisplacement(pI, disp);

    return disp;
}

void Foam::meshMetric::calculateGradients(const label& pI)
{
    const labelList& cl = mesh_.pointCells()[pI];
    const scalarField& AR = metric_.aspectRatio();
    derivative_ = vector::zero;
    hessian_ *= 0;

    calculateDerivatives(pI);
    forAll(cl, cI)
    {
           const scalar& volume = cellVols()[cl[cI]];
           const scalar& surface = cellSurface()[cl[cI]];
           const vector& d_volume = cellVolumeDerivatives()[cI];
           const vector& d_surface = cellSurfaceDerivatives()[cI];
           const tensor& d_2ndSurf = surfaceHessian()[cI];

           if (volume > ROOTVSMALL)
           {
            vector der = sphericityDer(volume, surface, d_volume, d_surface);
            tensor H = sphericity2ndDer(volume, surface, d_volume, d_surface, d_2ndSurf);
            der *= AR[cl[cI]]*AR[cl[cI]];
//            der /= cellTypeCoef()[cl[cI]];
            //modify hessian to minimize for "non-Sphericity"
            nonSphericityHessian(H, der, cI);
            hessian_ += H;
            derivative_ += -nonSphericity_[cI]/spher_[cI]*der;
           }
           else
           {
               derivative_ += d_volume;
           }
    }
    derivative_ /= pointCellsSize()[pI];
    hessian_ /= pointCellsSize()[pI];
}

void Foam::meshMetric::update(const label& pI)
{
    updateGeometry(pI);
}

void Foam::meshMetric::updateGlobal(const labelList& pList)
{
    forAll(pList, I)
    {
        const label& pI = pList[I];
        updateGeometry(pI);
    }
}

const Foam::label& Foam::meshMetric::currentLabel() const
{
    return pLabel_;
}

const Foam::vector& Foam::meshMetric::getDerivative(const label& pI) const
{
    if (pI != currentLabel())
    {
        Info<<"Asking for wrong derivative"<<nl
            <<"point inserted is not the current optimizing point"<<endl;
    }
    return derivative_;
}

const Foam::tensor& Foam::meshMetric::getHessian(const label& pI) const
{
    if (pI != pLabel_)
    {
        Info<<"Asking for wrong Hessian"<<nl
            <<"point inserted is not the current optimizing point"<<endl;
    }
    return hessian_;
}

void Foam::meshMetric::printInfo()
{
//    volScalarField Sph
//    (
//        IOobject
//        (
//            "Sph",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar("sph", dimless, 0)
//    );

    scalar minV = 2;
    scalar maxV = -1;
    forAll(mesh_.cells(), cI)
    {
        scalar sph = sphericityValue(cI);
/*
        if (sph>0.6*cellTypeCoef()[cI])
        {
            Info<<" Warning "<<nl
                <<" Cell Volume ="<<cellVols()[cI]<<nl
                <<" Cell Surface = "<<cellSurface()[cI]<<nl
                <<" Cell ID "<<cI<<endl;
        }
*/
        minV = min(sph, minV);
        maxV = max(sph, maxV);
    }

    reduce(
        std::tie(minV, maxV),
        ParallelOp<minOp<scalar>, maxOp<scalar>>{}
    );

    Info<<"Minimum sphericity value = "<<minV<<nl
        <<"Maximum sphericity value = "<<maxV<<endl;
}

// ************************************************************************* //
