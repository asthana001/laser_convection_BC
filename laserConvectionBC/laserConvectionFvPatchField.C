/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author
    Priyanshu Asthana
    Materials Modelling Group, IISc Bangalore
    India
    
\*---------------------------------------------------------------------------*/

#include "laserConvectionFvPatchField.H"
#include "scalarIOList.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::laserConvectionFvPatchField<Type>::laserConvectionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{
//    if (debug_)
    {
        Info<< "Constructor 1\n" << endl;
    }

    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::laserConvectionFvPatchField<Type>::laserConvectionFvPatchField
(
    const laserConvectionFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    nSources_(ptf.nSources_),
    reducedCoeff_(ptf.reducedCoeff_),
    reducedCoeffName_(ptf.reducedCoeffName_),
    coeff_(ptf.coeff_),
    kField_(ptf.kField_),
    kFieldName_(ptf.kFieldName_),
    kValue_(ptf.kValue_),
    sourceCenters_(ptf.sourceCenters_),
    normals_(ptf.normals_),
    sigmaX_(ptf.sigmaX_),
    sigmaY_(ptf.sigmaY_),
    power_(ptf.power_),
    correlation_(ptf.correlation_),
    HTCHeating_(ptf.HTCHeating_),
    HTCCooling_(ptf.HTCCooling_),
    TfluidHeating_(ptf.TfluidHeating_),
    TfluidQuench_(ptf.TfluidQuench_),
    heatingTime_(ptf.heatingTime_),
    heatingTimeSource_(ptf.heatingTimeSource_),
    motion_(ptf.motion_),
    motionMode_(ptf.motionMode_),
    startMotion_(ptf.startMotion_),
    centerUpdated_(ptf.centerUpdated_),
    actualTime_(ptf.actualTime_),
    motionCenters_(ptf.motionCenters_),
    startAngle_(ptf.startAngle_),
    omega_(ptf.omega_),
    nCycles_(ptf.nCycles_),
    lMPoints_(ptf.lMPoints_),
    timeAcc_(ptf.timeAcc_),
    endPoint_(ptf.endPoint_),
    gaussDistribution_(ptf.gaussDistribution_),
    heatFluxDistribution_(ptf.heatFluxDistribution_),
    valueFraction_(ptf.valueFraction_),
    refValue_(ptf.refValue_),
    debug_(ptf.debug_)
{
    if (debug_)
    {
        Info<< "Constructor 2\n" << endl;
    }
}


template<class Type>
Foam::laserConvectionFvPatchField<Type>::laserConvectionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    nSources_(0),
    reducedCoeff_(false),
    reducedCoeffName_("none"),
    coeff_(scalar(1)),
    kField_(false),
    kFieldName_("none"),
    kValue_(scalar(0)),
    sourceCenters_(0,point(0,0,0)),
    normals_(0,vector(0,0,0)),
    sigmaX_(0,scalar(0)),
    sigmaY_(0,scalar(0)),
    power_(0,scalar(0)),
    correlation_(0,scalar(0)),
    HTCHeating_(readScalar(dict.lookup("HTCheating"))),
    HTCCooling_(readScalar(dict.lookup("HTCquenching"))),
    TfluidHeating_(readScalar(dict.lookup("TfH"))),
    TfluidQuench_(readScalar(dict.lookup("TfQ"))),
    heatingTime_(readScalar(dict.lookup("heatingTime"))),
    heatingTimeSource_(0,scalar(0)),
    motion_(0,false),
    motionMode_(0,word("None")),
    startMotion_(0,scalar(0)),
    centerUpdated_(0,false),
    actualTime_(scalar(0.)),
    motionCenters_(0,point(0,0,0)),
    startAngle_(0,scalar(0)),
    omega_(0,scalar(0)),
    nCycles_(0,scalar(0)),
    lMPoints_(0,pointField(0,point(0,0,0))),
    lMSpeed_(0,scalar(0)),
    timeAcc_(0,scalarList(0,scalar(0))),
    endPoint_(0,false),
    gaussDistribution_(0,scalarField(0,scalar(0))),
    heatFluxDistribution_(p.size(), scalar(0)),
    valueFraction_(p.size(), scalar(1)),
    refValue_(TfluidHeating_),
    debug_(false)
{
    if (debug_)
    {
        Info<< "Constructor 3\n" << endl;
    }

    //- Name of the patch
    const word patchName = this->patch().name();

    //- TH::141116 Critical error based on non initiated field
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        FatalErrorIn
        (
            "laserConvectionFvPatchField.C\n"
            "value is missing in for patch " + patchName + "\n"
        ) << abort(FatalError);
    }

    Info<< "Reading laserConvection from patch " << patchName << "\n" << endl;

    //- First check if regIO of k and reduced coeff is available
    const wordList& regIO = this->db().names();

    //- Check if kFieldName is in the RegObj
    forAll(regIO, IO)
    {
        if (regIO[IO] == kFieldName_)
        {
            kField_ = true;

            Info<< " ++ k field " << kFieldName_ << " is used" << endl;
        }
    }

    //- Update
    reducedCoeffName_ = dict.lookupOrDefault<word>("powerReduceName", "none");
    kFieldName_ = dict.lookupOrDefault<word>("kValueName", "none");

    //- If kFieldName was not found in the RegObj and is not set to 'none'
    //  throw out an error message
    if (kField_ == false && kFieldName_ != "none")
    {
        forAll(regIO, IO)
        {
            Info<< "RegObj available: " << regIO[IO] << endl;
        }

        FatalErrorIn
        (
            "laserConvectionFvPatchField.C\n"
            "The kField that is given in boundary condition was not found.\n"
            + kFieldName_ + " is not in the registered objectives.\n"
            "You can use 'none' or a valid field name.\n"
        ) << abort(FatalError);
    }

    //- If no k field is specified, we need a value
    if (!kField_)
    {
        kValue_ = readScalar(dict.lookup("kValue"));

        Info<< " ++ No k field is used. Value for k is "
            << kValue_ << endl;
    }


    //- Bool variable if reducedCoeff is used
    //  false = not found not used
    //  true  = found and used
    bool rCNfound = false;

    forAll(regIO, IO)
    {
        //- Power-Reduce Coeff is used
        if (!rCNfound)
        {
            //- Explicit unset
            if (reducedCoeffName_ == "none")
            {
                //- Keyword found
                rCNfound = true;

                //- Leave the loop
                break;
            }
            else
            {
                //- Check if reduceCoeffName is in the database
                if (regIO[IO] == reducedCoeffName_)
                {
                    reducedCoeff_ = true;

                    rCNfound = true;

                    Info<< " ++ Power reduced coeff is used" << endl;

                    break;
                }
            }
        }
    }

    //- If no reducedCoeffName correlates give an ERROR
    if (!rCNfound)
    {
        //- Output the available fields
        forAll(regIO, IO)
        {
            Info<< regIO[IO] << "\n";
        }

        FatalErrorIn
        (
            "laserConvectionFvPatchField.C\n"
            "The reducedCoeffName field given in the boundary for patch "
            + patchName + " is not in the database.\n The"
            "valid fields are given above or you can use 'none'.\n"
            "Cannot find field '" + reducedCoeffName_ + "'\n"
        ) << abort(FatalError);
    }

    // Creating a dictionary "contentOfTable" with list of content
    const wordList& contentOfTable = dict.toc();

    //- Bool for source dict check; currently false but turns true if laser source is used.
    bool sourceDictFound = false;

    //- Iterating throught all the item of list declared above; 'c' : indexing variable 
    forAll(contentOfTable, c)
    {
        //- Checking if the current item is Subdict. If found -> new LASER source
        if (dict.isDict(contentOfTable[c])) //- ******* ERROR HERE ***** Not able to detect the source here.
        {
            //- New source found, number of sources increment
            nSources_++;

            //- Add new motion switch 
            motion_.append(false);

            //- Increase fields for motion and initialize
            motionMode_.append("off");
            startMotion_.append(0);
            motionCenters_.append(point(0,0,0));
            omega_.append(0);
            lMPoints_.append(pointField(0));
            lMSpeed_.append(0);
            timeAcc_.append(scalarField(0));
            endPoint_.append(false);
            centerUpdated_.append(false);
            startAngle_.append(0);
            nCycles_.append(0);

            //- Changing bool value (means at least 1 source exist)
            sourceDictFound = true;

            //- Adding name of source dict
            sourceDictName_.append(contentOfTable[c]);

            // Finds and returns the sub-dictionary, now sourceDict is a new dictionary
            const dictionary& sourceDict = dict.subDict(contentOfTable[c]);

            //- Add dictionary to field (for writing)
            sourceDicts_.append(sourceDict);

            //- Read center of beam
            sourceCenters_.append(sourceDict.lookup("spotCenter"));

            //- Read normal of beam
            normals_.append(sourceDict.lookup("normal"));

            //- Read some special stuff for the laser
            sigmaX_.append(readScalar(sourceDict.lookup("sigmaX")));
            sigmaY_.append(readScalar(sourceDict.lookup("sigmaY")));
            correlation_.append(readScalar(sourceDict.lookup("correlation")));
            power_.append(readScalar(sourceDict.lookup("power")));

            //- Special activation for the LASER source
            //  Thats for the following cases:
            //      The heat treatment is done with quenching but the laser
            //      should stop at t = t1 but the qenching procedure will
            //      start at t = t2. For that purpose, a heatingTimeSource_ is
            //      used to evaluate this problem.
            // Note: Not really needed but for completition
            const wordList& cOSD = sourceDict.toc();

            bool found = false;

            forAll(cOSD, w)
            {
                if (cOSD[w] == "active")
                {
                    found = true;

                    break;
                }
            }

            //- If active is found we set the LASER source time to this value
            if (found)
            {
                const bool active = sourceDict.lookup("active");

                if (active)
                {
                    heatingTimeSource_.append
                    (
                        readScalar(sourceDict.lookup("active"))
                    );
                }
            }
            //- Otherwise the LASER source will be active till the heating time is passed
            else
            {
                heatingTimeSource_.append(heatingTime_);
            }

            //- Motion section
            const wordList& contentOfSourceDict = sourceDict.toc();

            //- Simple read - no check if the dictionary is related to motion
            forAll(contentOfSourceDict, csd)
            {
                if (sourceDict.isDict(contentOfSourceDict[csd]))
                {
                    //- Get content of motion dictionary
                    const dictionary& motionDict =
                       sourceDict.subDict(contentOfSourceDict[csd]);

                    //- Set motion to true
                    motion_[nSources_-1] = true;

                    //- Start time for motion
                    startMotion_[nSources_-1] =
                        motionDict.lookupOrDefault("start", scalar(0));

                    //- Get the mode (linear or circular)
                    motionMode_[nSources_-1] = word(motionDict.lookup("mode"));

                    //- If motion mode is 'off', unset motion
                    if (motionMode_[nSources_-1] == "off")
                    {
                        motion_[nSources_-1] = false;
                    }
                }
            }

            //- If we found a subdictionary for the motion, read the stuff
            if (motion_[nSources_-1])
            {
                Info<< "\n ++ Read motion dict for source "
                    << sourceDictName_[nSources_-1] << endl;

                forAll(contentOfSourceDict, csd)
                {
                    if (sourceDict.isDict(contentOfSourceDict[csd]))
                    {
                        const dictionary& motionDict =
                           sourceDict.subDict(contentOfSourceDict[csd]);

                        //- Rotation motion
                        if (motionMode_[nSources_-1] == "circle")
                        {
                            Info<< "    ++ Circle mode for motion is used"
                                << endl;

                            motionCenters_[nSources_-1] =
                                motionDict.lookup("center");

                            nCycles_[nSources_-1] = motionDict.lookupOrDefault
                                (
                                    "nCycles",
                                    scalar(0)
                                );

                            omega_[nSources_-1] =
                                readScalar(motionDict.lookup("omega"));

                            endPoint_[nSources_-1] = false;

                            //- TH::151016 Take start angle into account

                            //- Vector between center of rotation and spot
                            const vector tmp =
                                sourceCenters_[nSources_-1]
                              - motionCenters_[nSources_-1];

                            //- Distance between center of rotation and spot
                            const scalar delta = mag(tmp);

                            //- Dx and dy
                            const scalar dxStart =
                                sourceCenters_[nSources_-1].x()
                              - motionCenters_[nSources_-1].x();

                            const scalar dyStart =
                                sourceCenters_[nSources_-1].y()
                              - motionCenters_[nSources_-1].y();

                            //- First and second quadrant dy positive
                            //  third and fourth quadrant dy negative

                            //- First and second
                            if (dyStart >= 0)
                            {
                                startAngle_[nSources_-1] =
                                    acos(dxStart / delta);
                            }

                            else
                            {
                                startAngle_[nSources_-1] =
                                  2*M_PI - acos(dxStart / delta);
                            }
                        }

                        //- Linear motion
                        else if (motionMode_[nSources_-1] == "linear")
                        {
                            Info<< "    ++ Linear mode for motion is used"
                                << endl;

                            lMSpeed_[nSources_-1] =
                               readScalar(motionDict.lookup("linearVelocity"));


                            if (motionDict.isDict("points"))
                            {
                                const dictionary& pDict =
                                    motionDict.subDict("points");

                                const int nPoints = pDict.toc().size();

                                //- Points have to be p(int) point coordinates
                                for (int p=0; p<nPoints; ++p)
                                {
                                    const word pointToEval =
                                        "p" + std::to_string(p);

                                    lMPoints_[nSources_-1].append
                                    (
                                        pDict.lookup(pointToEval)
                                    );
                                }
                            }

                            //- Set the spot center to the start point
                            sourceCenters_[nSources_-1] =
                                lMPoints_[nSources_-1][0];

                            //- Calculate the accumulation of time that is for the different linear parts | + length
                            for
                            (
                                int i=0;
                                i < lMPoints_[nSources_-1].size()-1;
                                ++i
                            )
                            {
                                //- Vector from p_n to p_(n+1)
                                const vector p1ToP2 =
                                    lMPoints_[nSources_-1][i+1]
                                  - lMPoints_[nSources_-1][i];

                                //- Length
                                const scalar length = mag(p1ToP2);

                                //- Time to reach the end point p_(n+1) from t = 0
                                if (i == 0)
                                {
                                    timeAcc_[nSources_-1].append
                                    (
                                        length / lMSpeed_[nSources_-1]
                                    );

                                    //lengthAcc_[nSources_-1].append(length);
                                }
                                else
                                {
                                    const scalar prevTime =
                                        timeAcc_[nSources_-1][i-1];

                                    timeAcc_[nSources_-1].append
                                    (
                                        prevTime
                                      +(length / lMSpeed_[nSources_-1])
                                    );
                                }
                            }
                        }
                        else if (motionMode_[nSources_-1] == "none")
                        {
                            Info<< "  ++ No motion mode is used \n" << endl;
                        }

                        //- Mode that is not available
                        else
                        {
                            FatalErrorIn
                            (
                                "laserConvectionFvPatchField.C\n"
                                "Mode motion is not defined. You try to use "
                                + motionMode_[nSources_-1] + "\nAvailable "
                                "modes are 'linear', 'circle' and 'none'"
                            ) << abort(FatalError);
                        }
                    }
                }
            }
        }
    }

    if (!sourceDictFound)
    {
        Info<< " ++ For patch " << patchName << " only convection"
            << " is applied, no laser source found.\n" << endl;
    }
    else
    {
        Info<< " ++ For patch " << patchName << " a laser source"
            << " is used. In total there are " << nSources_ << " sources\n"
            << endl;
    }

    //- Expand gaussDistribution_ field
    for(unsigned int i=0; i<nSources_; i++)
    {
        gaussDistribution_.append(scalarField(p.size()));
    }

    //- Build the refGrad, refValue and valueFraction for each face
    updateCoeffs();
}

template<class Type>
Foam::laserConvectionFvPatchField<Type>::laserConvectionFvPatchField
(
    const laserConvectionFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    nSources_(ptf.nSources_),
    reducedCoeff_(ptf.reducedCoeff_),
    reducedCoeffName_(ptf.reducedCoeffName_),
    coeff_(ptf.coeff_),
    kField_(ptf.kField_),
    kFieldName_(ptf.kFieldName_),
    kValue_(ptf.kValue_),
    sourceCenters_(ptf.sourceCenters_),
    normals_(ptf.normals_),
    sigmaX_(ptf.sigmaX_),
    sigmaY_(ptf.sigmaY_),
    power_(ptf.power_),
    correlation_(ptf.correlation_),
    HTCHeating_(ptf.HTCHeating_),
    HTCCooling_(ptf.HTCCooling_),
    TfluidHeating_(ptf.TfluidHeating_),
    TfluidQuench_(ptf.TfluidQuench_),
    heatingTime_(ptf.heatingTime_),
    heatingTimeSource_(ptf.heatingTimeSource_),
    motion_(ptf.motion_),
    motionMode_(ptf.motionMode_),
    startMotion_(ptf.startMotion_),
    centerUpdated_(ptf.centerUpdated_),
    actualTime_(ptf.actualTime_),
    motionCenters_(ptf.motionCenters_),
    startAngle_(ptf.startAngle_),
    omega_(ptf.omega_),
    nCycles_(ptf.nCycles_),
    lMPoints_(ptf.lMPoints_),
    timeAcc_(ptf.timeAcc_),
    endPoint_(ptf.endPoint_),
    gaussDistribution_(ptf.gaussDistribution_),
    heatFluxDistribution_(ptf.heatFluxDistribution_),
    valueFraction_(ptf.valueFraction_),
    refValue_(ptf.refValue_),
    debug_(ptf.debug_)
{
    if (debug_)
    {
        Info<< "Constructor 5\n" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::laserConvectionFvPatchField<Type>::updateCoeffs()
{
    if (debug_)
    {
        Info<< "UpdateCoeffs()\n" << endl;
    }

    if (this->updated())
    {
        return;
    }

    const scalar& time = this->db().time().value();

    //- Update time if we are in a new timestep
    if (actualTime_ != time)
    {
        actualTime_ = time;

        //- Reset bools
        forAll(centerUpdated_, i)
        {
            centerUpdated_[i] = false;
        }
    }

    //- Update spot center
    updateSpotCenter();

    //- Update the valueFraction
    updateValueFraction();

    //- Update the reference temperature (refValue_)
    updateRefTemperature();

    //- Check if each single laser source is finished
    bool allEndPointsReached = true;

    forAll(endPoint_, i)
    {
        if (!endPoint_[i])
        {
            allEndPointsReached = false;
            break;
        }
    }

    if (time > heatingTime_ || allEndPointsReached)
    {
        heatFluxDistribution_ = 0;
    }
    else
    {
        //- Update the gauss distribution
        updateGaussDistribution();

        //- Update the heatFlux distribution (refGrad_)
        updateHeatFluxDistribution();
    }

    //- Update
    if (debug_)
    {
        Info<< "Patch = " << this->patch().name()
            << "\n--------------------------------------------\n" << endl;

        //- Min/max refGrad
        const scalar minRefGrad = min(heatFluxDistribution_);
        const scalar maxRefGrad = max(heatFluxDistribution_);

        //- Min/max valueFraction
        const scalar minValueFraction = min(valueFraction_);
        const scalar maxValueFraction = max(valueFraction_);

        Info<< "Min refGrad = " << minRefGrad << endl;
        Info<< "Max refGrad = " << maxRefGrad << endl;

        Info<< "Min valueFraction = " << minValueFraction << endl;
        Info<< "Max valueFraction = " << maxValueFraction << endl;

        Info<< "refValue = " << refValue_ << endl;
    }

    this->refGrad() = pTraits<Type>::one * heatFluxDistribution_;
    this->refValue() = pTraits<Type>::one * refValue_;
    this->valueFraction() = valueFraction_;

    this->mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::laserConvectionFvPatchField<Type>::updateSpotCenter()
{
    const scalar& t = this->db().time().value();

    forAll(motion_, m)
    {
        if (debug_)
        {
            Info<< "Update spot center\n" << endl;
        }

        //- Motion only possible in x-y (z not considered)
        if (motion_[m] && !endPoint_[m])
        {
            //- If motion started and we need to update the center
            if (t > startMotion_[m] && !centerUpdated_[m])
            {
                //- Set center updated
                centerUpdated_[m] = true;

                //- Circular motion in x-y plane
                if (motionMode_[m] == "circle")
                {
                    const scalar angleRad = (t-startMotion_[m])*omega_[m];

                    //- Check if nCycles already done

                        //- Radian for nCyles 2pi = 1 cycle
                        const scalar radians = nCycles_[m] * 2 * M_PI;

                        //- Set power to zero and bool for end point reached
                        if (angleRad >= radians && radians > scalar(0))
                        {
                            endPoint_[m] = true;
                            power_[m] = 0;
                        }

                    //- Vector between two centers
                    const vector tmp = sourceCenters_[m] - motionCenters_[m];

                    //- Distance between two centers
                    const scalar delta = mag(tmp);

                    //- dx & dy
                    const scalar dx = cos(angleRad+startAngle_[m])*delta;
                    const scalar dy = sin(angleRad+startAngle_[m])*delta;

                    //- Update center coordinates
                    sourceCenters_[m][0] = motionCenters_[m].x()+dx;
                    sourceCenters_[m][1] = motionCenters_[m].y()+dy;
                }

                //- Linear motion in x-y plane
                else if (motionMode_[m] == "linear")
                {
                    //- Simple calculation of the laser motion. The
                    //  laser moves with the given speed, therefore
                    //  we calculate the vector inbetween the points
                    //  of id_ and id_+1. Out of the vector we build
                    //  the magnitude and check how long it takes to
                    //  move. After we reach the end point, the id_
                    //  is moved further.
                    label id = timeAcc_[m].size();

                    //- Find the points in which we are at the moment
                    forAll(timeAcc_[m], i)
                    {
                        if (timeAcc_[m][i] > t-startMotion_[m])
                        {
                            id = i;
                            break;
                        }
                    }

                    if (id >= timeAcc_[m].size())
                    {
                        endPoint_[m] = true;
                        power_[m] = 0;
                    }

                    //- Time acc | length acc
                    scalar timeReduce = 0;
                    //scalar lengthReduce = 0;

                    if (id != 0)
                    {
                        timeReduce = timeAcc_[m][id-1];
                    //    lengthReduce = lengthAcc_[m][id-1];
                    }

                    //- Vector from p1 to p2
                    const vector p1ToP2 =
                        lMPoints_[m][id+1] - lMPoints_[m][id];

                    //- Length
                    const scalar length = mag(p1ToP2);

                    //- Time to reach the end point p2
                    const scalar timeFromP1toP2 = length / lMSpeed_[m];

                    //- Ratio
                    const scalar ratio =
                        (t-timeReduce-startMotion_[m]) / timeFromP1toP2;

                    //- Update the center
                    sourceCenters_[m] = lMPoints_[m][id] + p1ToP2 * ratio;
                }
            }
        }
        else if (endPoint_[m])
        {
            Info<< " ++ laserConvectionFvPatchField.C\n"
                   " ++ There is no further linear point. End of "
                   "points reached. LASER turned off for "
                << sourceDictName_[m] << " on patch "
                << this->patch().name() << ".\n";
        }
    }
}


template<class Type>
void Foam::laserConvectionFvPatchField<Type>::updateValueFraction()
{
    if (debug_)
    {
        Info<< "Update valueFraction\n" << endl;
    }

    //- The current time
    const scalar& time = this->db().time().value();
    scalar HTC = 0;

    if (time > heatingTime_)
    {
        HTC = HTCCooling_;
    }
    else
    {
        HTC = HTCHeating_;
    }

    //- Inverse distance (face center - cell center)
    const scalarField& delta = this->patch().deltaCoeffs();

    scalarField kAtPatch(delta.size(), scalar(0));

    //- If k field is available, search and use
    if (kField_)
    {
        //- This patch
        const fvPatch& p = this->patch();

        //- Patch index
        const label& patchLabel = p.index();

        //- Thermal conductivity
        const volScalarField& kField =
            this->db().template lookupObject<volScalarField>(kFieldName_);

        kAtPatch = kField.boundaryField()[patchLabel];
    }
    else
    {
        kAtPatch = kValue_;
    }

    valueFraction_ = pow(kAtPatch*delta/HTC + 1, -1);
}


template<class Type>
void Foam::laserConvectionFvPatchField<Type>::updateHeatFluxDistribution()
{
    if (debug_)
    {
        Info<< "Update the heat flux distribution\n" << endl;
    }

    //- Set everything to zero
    heatFluxDistribution_ = 0;

    //- This patch
    const fvPatch& p = this->patch();

    scalarField kAtPatch(p.size(), scalar(0));

    //- If k field is available, search and use
    if (kField_)
    {

        //- Patch index
        const label& patchLabel = p.index();

        //- Thermal conductivity
        const volScalarField& kField =
            this->db().template lookupObject<volScalarField>(kFieldName_);

        kAtPatch = kField.boundaryField()[patchLabel];
    }
    else
    {
        kAtPatch = kValue_;
    }

    if (reducedCoeff_)
    {
        const scalarIOList& tmp_ =
            this->db().template lookupObject<scalarIOList>(reducedCoeffName_);

        coeff_ = tmp_[0];
    }

    const scalar& time = this->db().time().value();
    unsigned int active = 0;

    for(unsigned int i=0; i<nSources_; i++)
    {
        //- Check if source is active
        if (time > heatingTimeSource_[i])
        {
            active = 0;
        }
        else
        {
            active = 1;
        }

        heatFluxDistribution_ +=
           coeff_ * power_[i] * gaussDistribution_[i] / kAtPatch * active;
    }
}


template<class Type>
void Foam::laserConvectionFvPatchField<Type>::updateGaussDistribution()
{
    if (debug_)
    {
        Info<< "Update gauss distribution\n" << endl;
    }

    //- Patch
    const fvPatch& p = this->patch();

    //- Face centers
    const List<point>& cf = p.Cf();

    //- Face normals
    //const vectorField& fn = p.nf();

    for(unsigned int i=0; i<nSources_; i++)
    {
        //- Pre-Coefficient and pre-Coefficient for exponent
        const scalar pre =
            1/(2*M_PI*sigmaX_[i]*sigmaY_[i]
           *sqrt(1-sqr(correlation_[i])));
        const scalar expPre = -1/(2*(1-sqr(correlation_[i])));

        //  Laser only applicable on x-y at the moment
        //  TH::091116 Change names for better understanding
        //  Removed references
        const scalar xGaussCenter = sourceCenters_[i].x();
        const scalar yGaussCenter = sourceCenters_[i].y();

        //- Loop through all face centers
        forAll(cf, c)
        {
            //- Point coordinates of faces center
            const scalar xCf = cf[c][0];
            const scalar yCf = cf[c][1];

            const scalar X = xGaussCenter - xCf;
            const scalar Y = yGaussCenter - yCf;

            //- Pre-Coefficient
            const scalar expTerm =
                sqr(X)/sqr(sigmaX_[i])+sqr(Y)/sqr(sigmaY_[i])
               -(2*correlation_[i]*X*Y)/(sigmaX_[i]*sigmaY_[i]);

            gaussDistribution_[i][c] = pre * exp(expPre*expTerm);
        }
    }
}


template<class Type>
void Foam::laserConvectionFvPatchField<Type>::updateRefTemperature()
{
    if (debug_)
    {
        Info<< "Update reference temperature\n" << endl;
    }

    const scalar& time = this->db().time().value();

    if (time > heatingTime_)
    {
        refValue_ = TfluidQuench_;
    }
}


template<class Type>
void Foam::laserConvectionFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeKeyword("HTCheating");
    os<< HTCHeating_ << ";\n";

    os.writeKeyword("HTCquenching");
    os<< HTCCooling_ << ";\n";

    os.writeKeyword("TfH");
    os<< TfluidHeating_ << ";\n";

    os.writeKeyword("TfQ");
    os<< TfluidQuench_ << ";\n";
 
    os.writeKeyword("heatingTime");
    os<< heatingTime_ << ";\n";

    //- Write source dicts to next timestep
    forAll(sourceDicts_, dict)
    {
        os<< "\n";
        os.writeKeyword(sourceDictName_[dict]);
        os<< sourceDicts_[dict];
    }

    os.writeKeyword("powerReduceName");
    if (reducedCoeff_)
    {
        os<< reducedCoeffName_ << ";\n";
    }
    else
    {
        os<< "none;\n";
    }

    if (kField_)
    {
        os.writeKeyword("kName");
        os<< kFieldName_ << ";\n";
    }
    else
    {
        os.writeKeyword("kName none;\n");
        os.writeKeyword("kValue");
        os<< kValue_ << ";\n";
    }

    //- Write face values
    this->writeEntry("value", os); // old.OpenFOAM.org (<=6)
//    writeEntry(os, "value", *this); 	 // OpenFOAM.org (>=7)

}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::laserConvectionFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //
