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

#include "hessianMeshOptimization/meshOptimization/hessianMeshOptimization.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::hessianMeshOptimization::meshOptimization()
{
    Info<<"Initial Status "<<endl;
    metric_.printInfo();
    while (!methodPtr_->converge() && (currentIteration_ < numberOfIterations_))
    {
        Info<<"Iteration  "<<currentIteration_+1<<endl;
        label pointSize = methodPtr_->getActiveSet().activePoints().size();
        for (label pI=0; pI<pointSize; ++pI)
        {
            const label& cP = methodPtr_->getActiveSet().activePoints()[pI];

            if (!metric_.boundaryPoint()[cP])
            {
                optimizePoint(cP);
            }
        }

        if (returnReduce(metric_.processorPoints().size(), sumOp<label>()) != 0)
        {
            optimizeProcessorPoints();
        }
        currentIteration_++;
        metric_.printInfo();
    }
//            Time& runTime = const_cast<Time&>(mesh_.time());
//            runTime++;
//            mesh_.movePoints(state_);
//            mesh_.write();
    Info<<"Finished optimisation"<<endl;
}

void Foam::hessianMeshOptimization::optimizePoint(const label& cP)
{
    scalar objective = metric_.pointValue(cP);

    vector dS = metric_.displacement(cP);

    vector dSsmooth = smoother_.displacement(cP);

    if (currentIteration_ > smoothingIterationThreshold_)
    {
        satisfyWolfeConditions(cP, objective, dSsmooth);
    }
    else
    {
        satisfyWolfeConditions(cP, objective, dS);
    }
}

void Foam::hessianMeshOptimization::optimizeProcessorPoints()
{
    //Parallel handling (NEEDS REFACTORING)
    forAll(metric_.procPointsLists(), lI)
    {
        const labelList& spLabels = metric_.procPointsLists()[lI];
        scalarField objectiveValues(spLabels.size(), 0);
        vectorField derivativeList(spLabels.size(), vector::zero);
        tensorField hessianList(spLabels.size(), tensor::zero);
        scalarField maxDisp(spLabels.size(), 0);
        for (int I=0; I<spLabels.size(); I++)
        {
            const label& pI = spLabels[I];
            objectiveValues[I] = metric_.pointValue(pI);

            metric_.calculateGradients(pI);

            derivativeList[I] = metric_.getDerivative(pI);
            hessianList[I] = metric_.getHessian(pI);
        }
        metric_.processorPointsMaxDisplacement(maxDisp, spLabels);
        syncTools::syncPointList
        (
            mesh_,
            spLabels,
            objectiveValues,
            plusEqOp<scalar>(),
            scalar(0)
        );
        syncTools::syncPointList
        (
            mesh_,
            spLabels,
            derivativeList,
            plusEqOp<vector>(),
            vector::zero
        );
        syncTools::syncPointList
        (
            mesh_,
            spLabels,
            hessianList,
            plusEqOp<tensor>(),
            tensor::zero
        );
        vectorField disp(spLabels.size());
        for (int I=0; I<spLabels.size(); I++)
        {
            const vector& der = derivativeList[I];
            const tensor& H = hessianList[I];

            disp[I] = modifiedNewtonsMethod::newtonsStep(der, H);

            if (mag(disp[I])>maxDisp[I])
            {
                disp[I] *= maxDisp[I]/mag(disp[I]);
            }
            const label& pI = spLabels[I];
            if (metric_.boundaryProcessorPoint()[pI])
            {
                disp[I] = vector::zero;
            }
        }
        satisfyWolfeConditions
        (
            spLabels,
            objectiveValues,
            disp,
            derivativeList
        );
    }
}

bool Foam::hessianMeshOptimization::converge()
{
    //TODO convergence criterion
    return false;
}

void Foam::hessianMeshOptimization::satisfyWolfeConditions
(
    const label& pI,
    const scalar& initObjective,
    const vector& disp
)
{
    vector initialSt = state_[pI];
    state_[pI] -= disp;

    const vector& deriv = metric_.getDerivative(pI);

    bool wolfeConditions = false;
    scalar K0 = 0.;
    scalar K = 1;
    scalar coef = 0.0001;
    scalar directionDer = 0;
    scalar decrease = 0;
    label steps = 0;

    while (!wolfeConditions)
    {
        //Find the new objective with the modified K
        metric_.update(pI);
        scalar objective = metric_.pointValue(pI);
        wolfeConditions = true;
        //Directional Derivative
        directionDer = deriv&disp;
        scalar delta = coef*directionDer*K;
        //Store the old K value
        K0 = K;

        decrease = initObjective - mag(delta);
        bool sufficientDecrease = (objective<=decrease);
        if (!sufficientDecrease)
        {
            //Set K to have the value that minimizes the quadratic function
            //that passes from point pI with a value of objective and a derivative of
            //directionDerivative
            K = -(directionDer*sqr(K0))/(2*(objective-initObjective-K0*directionDer));
            //Flags to speed up the correction if the quadratic is not a
            //good approximation of the current function
            if
            (
                (K > K0)
             || (mag(K0)-mag(K)<0.001)
             || (K<0)
            )
            {
                K = 0.5*K0;
            }
            wolfeConditions = false;
        }

        if (!wolfeConditions)
        {
            if (steps > 4)
            {
                K = 0;
                wolfeConditions = true;
            }
            state_[pI] = initialSt - K*disp;
            steps++;
        }
    }
}

void Foam::hessianMeshOptimization::satisfyWolfeConditions
(
    const labelList& pList,
    const scalarList& initObjective,
    const vectorField& disp,
    const vectorField& derivatives
)
{
    vectorField initialSt(pList.size());

    for (int I=0; I<pList.size(); I++)
    {
        const label& pI = pList[I];
        initialSt[I] = state_[pI];
        state_[pI] -= disp[I];
    }

    bool wolfeConditions = false;
    List<bool> armijo(pList.size(), false);

    scalarField K0(pList.size(), 0);
    scalarField K(pList.size(), 1);
    scalar directionDer = 0;
    scalar coef = 0.0001;
    scalar decrease = 0;
    label steps = 0;

    while (!wolfeConditions)
    {
        //Find the new objective with the modified K
        metric_.updateGlobal(pList);
        scalarField objective = metric_.pointValueGlobal(pList)();
        wolfeConditions = true;

        forAll(pList, pI)
        {
            if (!armijo[pI])
            {
                //Directional Derivative
                directionDer = derivatives[pI]&disp[pI];
                scalar delta = coef*directionDer*K[pI];
                //Store the old K value
                K0[pI] = K[pI];

                decrease = initObjective[pI] - mag(delta);
                bool sufficientDecrease = (objective[pI]<=decrease);
                if (!sufficientDecrease)
                {
                    //Set K to have the value that minimizes the quadratic function
                    //that passes from point pI with a value of objective and a derivative of
                    //directionDerivative
                    K[pI] = -(directionDer*sqr(K0[pI]))
                        /(2*(objective[pI]-initObjective[pI]-K0[pI]*directionDer));
                    //Flags to speed up the correction if the quadratic is not a
                    //good approximation of the current function
                    if
                    (
                        (K[pI] > K0[pI])
                     || (mag(K0[pI])-mag(K[pI])<0.001)
                     || (K[pI]<0)
                    )
                    {
                        K[pI] = 0.5*K0[pI];
                    }
                    wolfeConditions = false;
                }
                else
                {
                    armijo[pI] = true;
                }
            }
        }
        reduce(wolfeConditions, minOp<bool>());

        if (!wolfeConditions)
        {
            if (steps > 4)
            {
                forAll(pList, pI)
                {
                    if (!armijo[pI])
                    {
                        K[pI] = 0;
                    }
                }
                wolfeConditions = true;

            }
            forAll(pList, pI)
            {
                const label& pL = pList[pI];
                state_[pL] = initialSt[pI] - K[pI]*disp[pI];
            }
            steps++;
        }
        reduce(wolfeConditions, minOp<bool>());
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hessianMeshOptimization::hessianMeshOptimization
(
    fvMesh& mesh,
    const pointField& state,
    const dictionary& dict,
    const bool local
)
:
    dict_(dict),
    mesh_(mesh),
    state_(state),
    metric_(mesh_, state_, dict_),
    methodPtr_
    (
        local ?
        meshOptimizationMethod::New(dict_, mesh_, metric_)
        : meshOptimizationMethod::NewFull(dict_, mesh_, metric_)
    ),
    numberOfIterations_(dict_.lookupOrDefault("iterations", 10)),
    smoother_(metric_),
    currentIteration_(0),
    smoothingIterationThreshold_
    (
        dict_.lookupOrDefault<label>("smoothingStartingIteration", 10e8)
    )
{
}

Foam::hessianMeshOptimization::hessianMeshOptimization
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool local
)
:
    dict_(dict),
    mesh_(mesh),
    state_(mesh.points()),
    metric_(mesh_, state_, dict_),
    methodPtr_
    (
        local ?
        meshOptimizationMethod::New(dict_, mesh_, metric_)
        : meshOptimizationMethod::NewFull(dict_, mesh_, metric_)
    ),
    numberOfIterations_(dict_.lookupOrDefault("iterations", 10)),
    smoother_(metric_),
    currentIteration_(0),
    smoothingIterationThreshold_
    (
        dict_.lookupOrDefault<label>("smoothingStartingIteration", 10e8)
    )
{
}

Foam::hessianMeshOptimization::hessianMeshOptimization
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    dict_(dict),
    mesh_(mesh),
    state_(mesh.points()),
    metric_(mesh_, state_, dict_),
    methodPtr_
    (
        meshOptimizationMethod::New(dict_, mesh_, metric_)
    ),
    numberOfIterations_(dict_.lookupOrDefault("iterations", 10)),
    smoother_(metric_),
    currentIteration_(0),
    smoothingIterationThreshold_
    (
        dict_.lookupOrDefault<label>("smoothingStartingIteration", 10e8)
    )
{
}

Foam::autoPtr<Foam::hessianMeshOptimization>
Foam::hessianMeshOptimization::fullMeshOptimization
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    autoPtr<hessianMeshOptimization> optimPtr
        (new hessianMeshOptimization(mesh, dict, true));
    return optimPtr;
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::hessianMeshOptimization::newPoints()
{
    meshOptimization();
    return state_;
}

Foam::tmp<Foam::pointField> Foam::hessianMeshOptimization::newPoints
(
    const labelList& pointSet
)
{
    methodPtr_->merge(pointSet);
    meshOptimization();
    tmp<pointField> tnewPoints
    (
        new pointField(pointSet.size(), vector::zero)
    );
    pointField& newPoints = tnewPoints.ref();
    forAll(newPoints, pI)
    {
        newPoints[pI] = state_[pointSet[pI]];
    }
    return tnewPoints;
}
// ************************************************************************* //
