/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::laserConvectionFvPatchField

Description
    This boundary condition provides a laser source boundary (multi laser
    sources possible) and a convective heat transfer with HTC. If no source
    is added, the boundary works as a convective bc.

    HTCheating          HTC coeff during heating [W/m^2/K]
    HTCquencing         HTC coeff during cooling [W/m^2/K]
    TfH                 Temperature of surouding fluid during heating [K]
    TfQ                 Temperature of fluid during Quenching [K]

    sourceXY { }        Definition of a Gaussian Laser source
    power               power [W]
    spotCenter          center of gaussian
    normal              direction of the gaussian
    sigmaX              streching in x [m]
    sigmaY              streching in y [m]
    correlation         modification to shape
    mu                  not-used

    heatingTime         time of laser source active

    arbitraryName{}     Each subdictionary in the boundary condition defines a
                        own laser source. Therefore you can add as much as you
                        like to have.

    motion { }          Definition of the laser motion
                        NOTE: Keep in mind that the BC is fully explicit, that
                        means that it might be necessary to use very small
                        time steps in order to get the correct result.

    mode                mode of laser motion
                        mode options: circle, linear, off
    center              center of circle motion
    omega               rotation speed in [rad/s]
    nCycles             cycles to run (0 -> do not stop)

    points { }          Definition of points for linear motion
    linearVelocity      linear speed for the laser for linear motion [m/s]

    active              activate or deactivate the laser source

    powerReduceName     name of field (scalar) that reduce the LASER power
                        If the solver offers a registered scalar that can be
                        looked up from the db(), you can modify this scalar
                        from 0-1 within your solver to adjust your LASER.
                        The scalar will be multiplied by the power of the
                        LASER. Therefore you can make a temperature depended
                        energy insertion.
                        (Either you compile the bc with your solver or you
                         make sure that you have the object in the database)

    kName               temperature depended thermal conductivity field
                        If not available set to none and use kValue

    kValue              constant thermal conductivity [W/mK]

    Examples of the boundary condition specification

    Simple convection
    \verbatim
    myPatch
    {
        type            laserConvection;
        value           uniform 300;

        HTCheating      150;
        HTCquenching    15000;

        TfH             300;
        TfQ             300;

        heatingTime     10;
        kValue          24;
    }


    \endverbatim

    \verbatim
    myPatch
    {
        type            laserConvection;
        value           uniform 300;

        HTCheating      23;
        HTCquenching    15000;

        TfH             300;
        TfQ             290;

        heatingTime     100;

        source1
        {
            power       1000;
            spotCenter  (0.1 0 0);
            normal      (0 1 0);  //- Not used

            sigmaX      0.03;
            sigmaY      0.03;

            correlation 0;

            motion
            {
                mode    circle;
                center  (0 0 0);
                omega   6;
                nCycles 1;
                active  true;

                start   0.3;

                //mode linear;
                //points
                //{
                //  p0  (0 0 0);
                //  p1  (0 1 0);
                //  p2  (1 2 0);
                //}
                //linearVelocity 0.4;
            }
        }

        sourcexY
        {
        }

        LaserForX-Motion
        {
        }

        powerReduceName   none;
        kName             k;

        //kValue          250;
    }
    \endverbatim

SeeAlso
    Foam::mixedFvPatchField

SourceFiles
    laserConvectionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef laserConvectionFvPatchField_H
#define laserConvectionFvPatchField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class laserConvectionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class laserConvectionFvPatchField
:
    public mixedFvPatchField<Type>
{

private:

    // Debug switch

    // Member data

        //- Number of sources
        unsigned int nSources_;

        //- Source dictionaries
        List<dictionary> sourceDicts_;

        //- Name of source dict
        List<word> sourceDictName_;

        //- Solver has reduced coeff
        bool reducedCoeff_;

        //- Name of reduced power coeff
        word reducedCoeffName_;

        //- Power reduced coeff
        scalar coeff_;

        //- Lambda field available
        bool kField_;

        //- Lambda field name
        word kFieldName_;

        //- Lambda value
        scalar kValue_;

        //- List of source centers
        List<point> sourceCenters_;

        //- List of normals
        List<vector> normals_;

        //- List of sigmaX and sigY [1/m^2]
        List<scalar> sigmaX_;
        List<scalar> sigmaY_;

        //- List of power [W]
        List<scalar> power_;

        //- List of correlation factors
        List<scalar> correlation_;

        //- HTC coefficents
        scalar HTCHeating_;
        scalar HTCCooling_;

        //- Temperature of sourronding liquid during heating and quenching
        scalar TfluidHeating_;
        scalar TfluidQuench_;

        //- Heating time
        scalar heatingTime_;

        //- Laser active time
        List<scalar> heatingTimeSource_;


    // Motion stuff

        //- Motion switch
        List<bool> motion_;

        //- Mode of motion
        List<word> motionMode_;

        //- Start motion
        List<scalar> startMotion_;

        //- For pimple mode (center updated)
        List<bool> centerUpdated_;

        //- Actual time
        scalar actualTime_;

        //- Rotation around a center

            //- Center of motion
            List<point> motionCenters_;

            //- Start angle
            List<scalar> startAngle_;

            //- Speed in RAD
            List<scalar> omega_;

            //- n cyclics
            List<scalar> nCycles_;

        //- Linear motion

            //- Point list for linear motion
            List<pointField> lMPoints_;

            //- Speed of linear motion [m/s]
            List<scalar> lMSpeed_;

            //- Time accumulated between the points [s]
            List<scalarList> timeAcc_;

        // //- Rotation in a spiral manner

        //     //- Coeff to make it spiral
        //     List<scalar> spiralcoeff_;

        //- Turn laser off after end point reached
        List<bool> endPoint_;


    // For access

        //- Gauss distribution of the laser sources [1/m^2]
        List<scalarField> gaussDistribution_;

        //- refGrad
        //- HeatFlux distribution of the laser sources [W/m^2]
        scalarField heatFluxDistribution_;

        //- valueFraction
        scalarField valueFraction_;

        //- refValue
        scalar refValue_;

        //- patch values
        scalarField patchValues_;


    // Debug

        bool debug_;


public:

    //- Runtime type information
    TypeName("laserConvection");


    // Constructors

        //- Construct from patch and internal field
        laserConvectionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        laserConvectionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given laserConvectionFvPatchField onto a new patch
        laserConvectionFvPatchField
        (
            const laserConvectionFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        laserConvectionFvPatchField
        (
            const laserConvectionFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new laserConvectionFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        laserConvectionFvPatchField
        (
            const laserConvectionFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new laserConvectionFvPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Update spot center
        virtual void updateSpotCenter();

        //- Update value fraction
        virtual void updateValueFraction();

        //- Update gauss distribution based on input file
        virtual void updateGaussDistribution();

        //- Update the heat flux distribution (refGrad)
        virtual void updateHeatFluxDistribution();

        //- Update the reference temperature (refValue)
        virtual void updateRefTemperature();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<Type>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "laserConvectionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
