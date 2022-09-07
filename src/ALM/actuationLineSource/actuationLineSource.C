/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "actuationLineSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuationLineSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuationLineSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::actuationLineSource::printDebug()
{
    Info<< "Print Debugging Information" << endl;
    Info<< "turbineType = " << turbineType << endl;
    Info<< "baseLocation = " << baseLocation << endl;
    Info<< "numBladePoints = " << numBladePoints << endl;
    Info<< "pointDistType = " << pointDistType << endl;
    Info<< "sampleDistScalar = " << sampleDistScalar << endl;
    Info<< "sampleDistParam = " << sampleDistParam << endl;
    Info<< "epsilon = " << epsilon << endl;
    Info<< "smearRadius = " << smearRadius << endl;
    Info<< "sphereRadiusScalar = " << sphereRadiusScalar << endl;
    Info<< "azimuth = " << azimuth << endl;
    Info<< "time = " << runTime_.value() << endl;
    Info<< "pitch = " << pitch << endl;
    Info<< "rotSpeed = " << rotSpeed << endl;
    Info<< "nacYaw = " << nacYaw << endl << endl << endl;
    Info<< "numTurbinesDistinct = " << numTurbinesDistinct << endl;
    Info<< "turbineTypeDistinct = " << turbineTypeDistinct << endl;
    Info<< "turbineTypeID = " << turbineTypeID << endl << endl << endl;;
    Info<< "NumBl = " << NumBl << endl;
    Info<< "TipRad = " << TipRad << endl;
    Info<< "HubRad = " << HubRad << endl;
    Info<< "UndSling = " << UndSling << endl;
    Info<< "OverHang = " << OverHang << endl;
    Info<< "TowerHt = " << TowerHt << endl;
    Info<< "Twr2Shft = " << Twr2Shft << endl;
    Info<< "ShftTilt = " << ShftTilt << endl;
    Info<< "PreCone = " << PreCone << endl;
    Info<< "PitchRate = " << PitchRate << endl;
    Info<< "PitchAmplitude = " << PitchAmplitude << endl;
    Info<< "PitchControllerType = " << PitchControllerType << endl;
    Info<< "AirfoilType = " << AirfoilType << endl;
    Info<< "BladeData = " << BladeData << endl;
    Info<< "BladeStation = " << BladeStation << endl;
    Info<< "BladeP1 = " << BladeP1 << endl;
    Info<< "BladeP2 = " << BladeP2 << endl;
    Info<< "AirfoilTypesDistinct = " << airfoilTypesDistinct << endl;
    Info<< "BladeAirfoilTypeID = " << BladeAirfoilTypeID << endl << endl
        << endl;
    Info<< "airfoilAlpha = " << airfoilAlpha << endl;
    Info<< "airfoilCl = " << airfoilCl << endl;
    Info<< "airfoilCd = " << airfoilCd << endl;
    Info<< "db = " << db << endl;
    Info<< "bladePoints = " << bladePoints << endl;
    Info<< "bladeRadius = " << bladeRadius << endl;
    Info<< "towerShaftIntersect = " << towerShaftIntersect << endl;
    Info<< "rotorApex = " << rotorApex << endl;
    Info<< "uvShaft = " << uvShaft << endl;
    Info<< "uvShaftDir = " << uvShaftDir << endl;
    Info<< "uvTower = " << uvTower << endl;
    Info<< "deltaNacYaw = " << deltaNacYaw << endl;
    Info<< "deltaAzimuth = " << deltaAzimuth << endl;
    Info<< "bladeAlignedVectors = " << bladeAlignedVectors << endl;
    Info<< "chordAlignedVectors = " << chordAlignedVectors << endl;
    Info<< "samplePoints = " << samplePoints << endl;
    Info<< "bladeForce = " << bladeForce << endl;
}


void Foam::fv::actuationLineSource::computeRotSpeed()
{
    // Proceed turbine by turbine.
    forAll(rotSpeed, i)
    {
        deltaAzimuth[i] = rotSpeed[i] * dt;
    }
}


void Foam::fv::actuationLineSource::computePitch()
{
    // Proceed turbine by turbine.
    forAll(pitchZero, i)
    {
        int j = turbineTypeID[i];
        if (PitchControllerType[j] == "none")
            {
                pitch[i] = pitchZero[i];
            }
        else if (PitchControllerType[j] == "sinusoidal")
            {
                pitch[i] =
                    pitchZero[i] + PitchAmplitude[i]
                   *sin(PitchRate[i]*runTime_.value());
            }
    }
}



void Foam::fv::actuationLineSource::rotateBlades()
{
    // Perform rotation turbine by turbine.
    forAll(uvShaft, i)
    {
        // Check the rotation direction first and set the local delta azimuth
        // variable accordingly.
        scalar deltaAzimuthI = 0.0;
        if(restart_flag[i] == 1)
        {
            if (rotationDir[i] == "cw")
            {
                deltaAzimuthI =  azimuth[i] + rotSpeed[i] * runTime_.value();
            }
            if (rotationDir[i] == "ccw")
            {
                deltaAzimuthI = -azimuth[i] - rotSpeed[i] * runTime_.value();
            }
            restart_flag[i] = 0;
        }
        else
        {
            if (rotationDir[i] == "cw")
            {
                deltaAzimuthI =  deltaAzimuth[i];
            }
            if (rotationDir[i] == "ccw")
            {
                deltaAzimuthI = -deltaAzimuth[i];
            }
        }

        // Rotate turbine blades, blade by blade, point by point.
        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                bladePoints[i][j][k] =
                    rotatePoint
                    (
                        bladePoints[i][j][k],
                        rotorApex[i],
                        uvShaft[i],
                        deltaAzimuthI
                    );
            }
        }

        // Calculate the new azimuth angle and make sure it isn't bigger than
        // 2*pi.
        azimuth[i] = azimuth[i] + deltaAzimuthI;
        if (azimuth[i] > 2.0 * Foam::constant::mathematical::pi)
        {
            azimuth[i] -= 2.0 * Foam::constant::mathematical::pi;
        }
    }
}


void Foam::fv::actuationLineSource::computeNacYaw()
{
    // Proceed turbine by turbine.
    forAll(deltaNacYaw, i)
    {
        deltaNacYaw[i] = 0.0;
    }
}


void Foam::fv::actuationLineSource::yawNacelle()
{
    // Perform rotation turbine by turbine.
    forAll(uvTower, i)
    {
        // Rotate the rotor apex first.
        rotorApex[i] =
            rotatePoint
            (
                rotorApex[i],
                towerShaftIntersect[i],
                uvTower[i],
                deltaNacYaw[i]
            );

        // Recompute the shaft unit vector since the shaft has rotated.
        uvShaft[i] = rotorApex[i] - towerShaftIntersect[i];
        uvShaft[i] = (uvShaft[i]/mag(uvShaft[i])) * uvShaftDir[i];

        // Rotate turbine blades, blade by blade, point by point.
        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                bladePoints[i][j][k] =
                    rotatePoint
                    (
                        bladePoints[i][j][k],
                        towerShaftIntersect[i],
                        uvTower[i],
                        deltaNacYaw[i]
                    );
            }
        }

        // Compute the new yaw angle and make sure it isn't bigger than 2*pi.
        nacYaw[i] = nacYaw[i] + deltaNacYaw[i];
        if (nacYaw[i] >= 2.0 * Foam::constant::mathematical::pi)
        {
            nacYaw[i] -= 2.0 * Foam::constant::mathematical::pi;
        }
    }
}


void Foam::fv::actuationLineSource::computeChordAlignedVectors()
{
    for (int i = 0; i < numTurbines; i++)
    {
        int m = turbineTypeID[i];
        for (int j = 0; j < NumBl[m]; j++)
        {
            // If clockwise rotating, this vector points along the blade toward
            // the tip. If counter-clockwise rotating, this vector points along
            // the blade toward the root.
            if (rotationDir[i] == "cw")
            {
                bladeAlignedVectors[i][j][2] =
                    bladePoints[i][j][0] - rotorApex[i];
                bladeAlignedVectors[i][j][2] =
                    bladeAlignedVectors[i][j][2]
                   /mag(bladeAlignedVectors[i][j][2]);
            }
            else if (rotationDir[i] == "ccw")
            {
                bladeAlignedVectors[i][j][2] =
                    -(bladePoints[i][j][0] - rotorApex[i]);
                bladeAlignedVectors[i][j][2] =
                    bladeAlignedVectors[i][j][2]
                   /mag(bladeAlignedVectors[i][j][2]);
            }

            // This vector points in the tangential direction opposite the
            // turbines rotation type.  It is set up this way because it will
            // point in the direction of oncoming flow that the blade sees due
            // to rotation.
            bladeAlignedVectors[i][j][1] =
                bladeAlignedVectors[i][j][2] ^ uvShaft[i];
            bladeAlignedVectors[i][j][1] =
                bladeAlignedVectors[i][j][1]/mag(bladeAlignedVectors[i][j][1]);

            // This vector points normal to the other two and toward downwind
            // (not exactly downwind if the blade is coned).  It points in the
            // direction of the oncoming flow due to wind that the blade sees.
            bladeAlignedVectors[i][j][0] =
                bladeAlignedVectors[i][j][1]^bladeAlignedVectors[i][j][2];
            bladeAlignedVectors[i][j][0] =
                bladeAlignedVectors[i][j][0]/mag(bladeAlignedVectors[i][j][0]);

            for (int k = 0; k < numBladePoints[i]; k++)
            {
                // Rotate bladeAlignedVectors about (their) z axis to give
                // chordAlignedVectors
                vector axis = bladeAlignedVectors[i][j][2];
                // clockwise rotation about z-bladeAlignedVector
                scalar angle = -(Twist[i][j][k] + pitch[i])*degRad;
                for (int q = 0; q < 3; q++)
                {
                    chordAlignedVectors[i][j][k][q] =
                        rotatePoint
                        (
                            bladeAlignedVectors[i][j][q],
                            vector::zero,
                            axis,
                            angle
                        );
                }
            }
        }
    }
}


void Foam::fv::actuationLineSource::computeSamplePoints()
{
    for (int i = 0; i < numTurbines; i++)
    {
        int m = turbineTypeID[i];
        for (int j = 0; j < NumBl[m]; j++)
        {
            for (int k = 0; k < numBladePoints[i]; k++)
            {
                if (sampleDistParam[i] >= 0.0)
                {
                    rs[i][j][k] = sampleDistParam[i];
                }
                else
                {
                    rs[i][j][k] = sampleDistScalar[i]*smearingRadiusC[i][j][k];
                }

                // Sample points build-up
                if(samplingMethod[i] == "Schluntz3pointConstZ")
                {
                    // sample points at gamma = pi/2, pi and 3/2pi, where gamma
                    // is the azimuth around the blade starting from the trailing
                    // edge. Sample points at 1 chord from actuator line.
                     samplePoints[i][j][k][0] =
                         bladePoints[i][j][k]
                       + chordAlignedVectors[i][j][k][0]*rs[i][j][k];
                     samplePoints[i][j][k][1] =
                         bladePoints[i][j][k]
                       - chordAlignedVectors[i][j][k][1]*rs[i][j][k];
                     samplePoints[i][j][k][2] =
                         bladePoints[i][j][k]
                       - chordAlignedVectors[i][j][k][0]*rs[i][j][k];
                }
                else if(samplingMethod[i] == "Schluntz3point")                    
                {
                    // sample points at gamma = pi/2, pi and 3/2pi, where gamma
                    // is the azimuth around the blade starting from the trailing
                    // edge. Sample points at 1 chord from actuator line.
                    vector  ds = vector::zero;
                    scalar   r=0.0;
                    scalar  dl=0.0;
                    scalar phi=0.0;

                    for (int iP = 0; iP < numSamplePoints[i]; iP++)
                    {
                        if(iP == 0)
                        {
                            ds = chordAlignedVectors[i][j][k][0]*rs[i][j][k];
                        }
                        else if(iP == 1)
                        {
                            ds = -chordAlignedVectors[i][j][k][1]*rs[i][j][k];
                        }
                        else if(iP == 2)
                        {
                            ds = - chordAlignedVectors[i][j][k][0]*rs[i][j][k];
                        }
                        
                          r = sqrt(pow(bladePoints[i][j][k][1],2.0) + pow(bladePoints[i][j][k][2],2.0));
                         dl = sqrt(pow(ds[1],2.0) + pow(ds[2],2.0));
                        phi = dl/r;
                       
                        // decide the sign of phi, cw: positive, ccw: negative
                        if ( fabs(bladePoints[i][j][k][1]) >= fabs(bladePoints[i][j][k][2]) )
                        {
                            if (bladePoints[i][j][k][1] * ds[2] > 0.0)
                            {
                                phi = 1.0*phi;
                            }
                            else
                            {
                                phi = -1.0*phi;
                            }
                        }
                        else
                        {
                            if (bladePoints[i][j][k][2] * ds[1] < 0.0)
                            {
                                phi = 1.0*phi;
                            }
                            else
                            {
                                phi = -1.0*phi;
                            }
                        }

                        // rotate blade points to get sample points at the same radial position
                        samplePoints[i][j][k][iP] =
                        rotatePoint
                        (
                            bladePoints[i][j][k],
                            rotorApex[i],
                            uvShaft[i],
                            phi
                        );
                        // apply streamwise shift to the sample points
                        samplePoints[i][j][k][iP][0] = bladePoints[i][j][k][0] + ds[0];
                    }
                }
                else if (samplingMethod[i] == "lineAverage")
                {
                    // sample points at arbitrary equal sections of the chord
                    vector  ds = vector::zero;
                    scalar   r=0.0;
                    scalar  dl=0.0;
                    scalar phi=0.0;
                    scalar rotAng = 2*Foam::constant::mathematical::pi/numSamplePoints[i];
                        
                    for (int iP = 0; iP < numSamplePoints[i]; iP++)
                    {
                        ds = rs[i][j][k]*
                            rotatePoint
                            (
                                -chordAlignedVectors[i][j][k][1],
                                vector::zero,
                                chordAlignedVectors[i][j][k][2],
                                iP*rotAng
                            );

                          r = sqrt(pow(bladePoints[i][j][k][1],2.0) + pow(bladePoints[i][j][k][2],2.0));
                         dl = sqrt(pow(ds[1],2.0) + pow(ds[2],2.0));
                        phi = dl/r;
                       
                        // decide the sign of phi, cw: positive, ccw: negative
                        if ( fabs(bladePoints[i][j][k][1]) >= fabs(bladePoints[i][j][k][2]) )
                        {
                            if (bladePoints[i][j][k][1] * ds[2] > 0.0)
                            {
                                phi = 1.0*phi;
                            }
                            else
                            {
                                phi = -1.0*phi;
                            }
                        }
                        else
                        {
                            if (bladePoints[i][j][k][2] * ds[1] < 0.0)
                            {
                                phi = 1.0*phi;
                            }
                            else
                            {
                                phi = -1.0*phi;
                            }
                        }

                        // rotate blade points to get sample points at the same radial position
                        samplePoints[i][j][k][iP] =
                        rotatePoint
                        (
                            bladePoints[i][j][k],
                            rotorApex[i],
                            uvShaft[i],
                            phi
                        );
                        // apply streamwise shift to the sample points
                        samplePoints[i][j][k][iP][0] = bladePoints[i][j][k][0] + ds[0];
                    }
                }
				else
				{
		            FatalErrorInFunction
		                << " Invalid sampling method "<<samplingMethod[i]<< "." 
		                << " Valid sampling methods are \n" 
		                << " lineAverage \n" 
		                << " Schluntz3point \n"  
		                << " Schluntz3pointConstZ \n"
		                << exit(FatalError);
				}
            }
        }
    }
}


void Foam::fv::actuationLineSource::computeBladeForce()
{
    samplePointsUChord = samplePointsU; // initialize size.
    forAll(samplePointsU, i)
    {
        int m = turbineTypeID[i];
        forAll(samplePointsU[i], j)
        {
            forAll(samplePointsU[i][j], k)
            {
                // Declare and define the rotation matrix.
                tensor RM;
                RM.xx() = chordAlignedVectors[i][j][k][0].x();
                RM.xy() = chordAlignedVectors[i][j][k][0].y();
                RM.xz() = chordAlignedVectors[i][j][k][0].z();
                RM.yx() = chordAlignedVectors[i][j][k][1].x();
                RM.yy() = chordAlignedVectors[i][j][k][1].y();
                RM.yz() = chordAlignedVectors[i][j][k][1].z();
                RM.zx() = chordAlignedVectors[i][j][k][2].x();
                RM.zy() = chordAlignedVectors[i][j][k][2].y();
                RM.zz() = chordAlignedVectors[i][j][k][2].z();

                // Angular velocity
                vector VR =
                    (rotSpeed[i]*bladeRadius[i][j][k]*cos(PreCone[m][j]))
                   *bladeAlignedVectors[i][j][1];

                forAll(samplePointsU[i][j][k], l) // three sample points
                {
                    // approximate as sample points are not at bladeRadius !
                    samplePointsU[i][j][k][l] += VR;
                    // Perform the rotation.
                    samplePointsUChord[i][j][k][l] =
                        RM & samplePointsU[i][j][k][l];
                }
            }
        }
    }

    forAll(qC, i)
    {
        int m = turbineTypeID[i];
        // Set the total thrust/torque of the turbine to zero.  Thrust/torque
        // will be summed on a blade-element-wise basis.
        thrust[i] = 0.0;
        torqueRotor[i] = 0.0;
        forAll(qC[i], j)
        {
            scalar p_gamma = 0.0;

            forAll(qC[i][j], k)
            {

                if ( samplingMethod[i] == "Schluntz3pointConstZ" || samplingMethod[i] == "Schluntz3point" )
                {
                    // Deduction of U_REL and alpha using circulation concept
                    qC[i][j][k] =
                        0.5
                       *(
                            samplePointsUChord[i][j][k][0].y()
                          - samplePointsUChord[i][j][k][2].y()
                        );
    
                    GammaC[i][j][k] =
                        2*Foam::constant::mathematical::pi*rs[i][j][k]*qC[i][j][k];

                    VmagC[i][j][k] =
                        Foam::sqrt
                        (
                            Foam::sqr
                            (
                                samplePointsUChord[i][j][k][1].x() - qC[i][j][k]
                            )
                          + Foam::sqr
                            (
                                samplePointsUChord[i][j][k][2].y() + qC[i][j][k]
                            )
                        );
                    alphaC[i][j][k] =
                        Foam::asin
                        (
                            (
                                samplePointsUChord[i][j][k][1].x() - qC[i][j][k]
                            )
                           /mag(VmagC[i][j][k])
                        )
                       /degRad;
                }
                else if ( samplingMethod[i] == "lineAverage" )
                {
                    vector samplePointsUChordAverage = vector::zero;
                    forAll(samplePointsUChord[i][j][k],iP)
                    {
                        samplePointsUChordAverage = samplePointsUChord[i][j][k][iP] + samplePointsUChordAverage;    
                    }
                    samplePointsUChordAverage = samplePointsUChordAverage/numSamplePoints[i];
                    VmagC[i][j][k]            = Foam::sqrt(Foam::sqr(samplePointsUChordAverage.x())+Foam::sqr(samplePointsUChordAverage.y()));
                    alphaC[i][j][k]           = Foam::asin(samplePointsUChordAverage.x()/VmagC[i][j][k])/degRad;

                    qC[i][j][k] = 0.0;
                    forAll(samplePointsUChord[i][j][k],iP)
                    {
                        qC[i][j][k] = qC[i][j][k] + fabs(mag(samplePointsUChord[i][j][k][iP]) - mag(samplePointsUChordAverage));
                    }
                    qC[i][j][k] = qC[i][j][k]/numSamplePoints[i];
    
                    GammaC[i][j][k] =
                        2*Foam::constant::mathematical::pi*rs[i][j][k]*qC[i][j][k];
                }
				else
				{
		            FatalErrorInFunction
		                << " Invalid sampling method "<<samplingMethod[i]<< "." 
		                << " Valid sampling methods are \n" 
		                << " lineAverage \n" 
		                << " Schluntz3point \n"  
		                << " Schluntz3pointConstZ \n"
		                << exit(FatalError);
				}

//              Calculation of lift and drag polars. 
//              Info << airfoilInterpolationMethod[m] << endl;
//              Method 1: Linear nterpolation Between lift and drag polars.
                if (airfoilInterpolationMethod[m] == "LinearInterpolation" )
                {
                    List<label> airfoilsMP = 
                        extrapolate
                        (
                            bladeRadius[i][j][k],
                            BladeStation[m],
                            BladeAirfoilTypeID[m]
                        );

                    DynamicList<label> BladeStationID(BladeStation.size());
                    for (int iID = 0; iID < BladeStation[m].size(); iID++)    
                    {
                        label a(iID);
                        BladeStationID.append(a);
                    }
                    List<label> BladeStationIDMP = 
                        extrapolate
                        (
                            bladeRadius[i][j][k],
                            BladeStation[m],
                            BladeStationID
                        );
                    DynamicList<scalar> clTemp(2);
                    DynamicList<scalar> cdTemp(2);
                    DynamicList<scalar> BladeStationTemp(2);
                    for ( int iTemp = 0; iTemp < 2; iTemp++ )
                    {
                        clTemp[iTemp] = 
                            interpolate
                            (
                                alphaC[i][j][k],
                                airfoilAlpha[airfoilsMP[iTemp]],
                                airfoilCl[airfoilsMP[iTemp]]
                            );
                        cdTemp[iTemp] = 
                            interpolate
                            (
                                alphaC[i][j][k],
                                airfoilAlpha[airfoilsMP[iTemp]],
                                airfoilCd[airfoilsMP[iTemp]]
                            );
                        BladeStationTemp[iTemp] = BladeStation[m][BladeStationIDMP[iTemp]];
                    }

                    if (BladeStationTemp[0] == BladeStationTemp[1])
                    {
                        Cl[i][j][k] = clTemp[0];
                        Cd[i][j][k] = cdTemp[0];
                    }
                    else
                    {
                        Cl[i][j][k] = 
                            interpolate
                            (    
                                bladeRadius[i][j][k],
                                BladeStationTemp,
                                clTemp
                            );

                        Cd[i][j][k] = 
                            interpolate
                            (    
                                bladeRadius[i][j][k],
                                BladeStationTemp,
                                cdTemp
                            );
                    }
//                  Info<<"LinearInterpolation: Colloc="<<k<<" r = "<<bladeRadius[i][j][k]<<" belongs in ["<<BladeStationTemp[0]<<","<<BladeStationTemp[1];
//                  Info<<"] , airfoil belongs in ["<<airfoilsMP[0]<<","<<airfoilsMP[1]<<"]"<<endl;
                }
                else if (airfoilInterpolationMethod[m] == "IDInterpolation" )
                {
                    // Find the local airfoil type.
                    label airfoil =
                        interpolate
                        (
                            bladeRadius[i][j][k],
                            BladeStation[m],
                            BladeAirfoilTypeID[m]
                        );

                    // Use airfoil look-up tables to get coefficient of lift and
                    // drag.
                    Cl[i][j][k] =
                        interpolate
                        (
                            alphaC[i][j][k],
                            airfoilAlpha[airfoil],
                            airfoilCl[airfoil]
                        );
                    Cd[i][j][k] =
                        interpolate
                        (
                            alphaC[i][j][k],
                            airfoilAlpha[airfoil],
                            airfoilCd[airfoil]
                        );
//            Info<<"IDInterpolation: Colloc="<<k<<" r = "<<bladeRadius[i][j][k];
//            Info<<"] , airfoil is closer to "<<airfoil<<endl;
                }
                else if (airfoilInterpolationMethod[m] == "TuInterpolation" )
                {
                    scalar Tu = rotorUturb[i]/mag(VmagC[i][j][k]);
//                  Info << Tu << endl;
                    List<label> airfoilsMP = 
                        extrapolate
                        (
                            Tu,
                            BladeP2[m],
                            BladeAirfoilTypeID[m]
                        );
//                  if(i == 0 && j == 0)
//                  {
//                      Info << bladeRadius[i][j][k] << "\t" << Tu << "\t" << airfoilsMP  << endl;
//                  }
                    DynamicList<label> BladeP2ID(BladeP2.size());
                    for (int iID = 0; iID < BladeP2[m].size(); iID++)    
                    {
                        label a(iID);
                        BladeP2ID.append(a);
                    }
//                  Info << "BladeP2ID = " << BladeP2ID << endl;
                    List<label> BladeP2IDMP = 
                        extrapolate
                        (
                            Tu,
                            BladeP2[m],
                            BladeP2ID
                        );
//                  Info << "BladeP2IDMP = " << BladeP2IDMP << endl;
                    DynamicList<scalar> clTemp(2);
                    DynamicList<scalar> cdTemp(2);
                    DynamicList<scalar> BladeP2Temp(2);
                    for ( int iTemp = 0; iTemp < 2; iTemp++ )
                    {
                        clTemp[iTemp] = 
                            interpolate
                            (
                                alphaC[i][j][k],
                                airfoilAlpha[airfoilsMP[iTemp]],
                                airfoilCl[airfoilsMP[iTemp]]
                            );
                        cdTemp[iTemp] = 
                            interpolate
                            (
                                alphaC[i][j][k],
                                airfoilAlpha[airfoilsMP[iTemp]],
                                airfoilCd[airfoilsMP[iTemp]]
                            );
                        BladeP2Temp[iTemp] = BladeP2[m][BladeP2IDMP[iTemp]];
                    }
//                  Info << "Tu = " << Tu << endl;
//                  Info << "BladeP2Temp \t clTemp \t cdTemp " << endl;
//                  for ( int iTemp = 0; iTemp < 2; iTemp++ )
//                  {
//                      Info << BladeP2Temp[iTemp] << "\t" << clTemp[iTemp] << "\t" << cdTemp[iTemp] << endl;
//                  }
                    if (BladeP2Temp[0] == BladeP2Temp[1])
                    {
                        Cl[i][j][k] = clTemp[0];
                        Cd[i][j][k] = cdTemp[0];
                    }
                    else
                    {
                        Cl[i][j][k] = 
                            interpolate
                            (    
                                Tu,
                                BladeP2Temp,
                                clTemp
                            );
                        Cd[i][j][k] = 
                            interpolate
                            (    
                                Tu,
                                BladeP2Temp,
                                cdTemp
                            );
                    }
//                  Info << "Cl = " << Cl[i][j][k] << endl;
//                  Info << "Cd = " << Cd[i][j][k] << endl;
                }
                else if (airfoilInterpolationMethod[m] == "TCInterpolation" )
                {
                    scalar TC = TCratio[i][j][k];
//                 scalar TC = 
//                     interpolate
//                     (
//                                 bladeRadius[i][j][k],
//                                 BladeStation[m],
//                                 BladeP2[m]
//                     );
//                 if(i == 0 && j == 0)
//                 {
//                 Info << bladeRadius[i][j][k] << "\t" << TC << endl;
//                 Info << "radius = " << bladeRadius[i][j][k] << endl;
//                 Info << "TC = " << TC << endl;
//                 }
                    List<label> airfoilsMP = 
                        extrapolate
                        (
                                    TC,
                                    BladeP2[m],
                                    BladeAirfoilTypeID[m]
                        );
//                 if(i == 0 && j == 0)
//                 {
//                 Info << bladeRadius[i][j][k] << "\t" << TC << "\t" << airfoilsMP  << endl;
//                 }
                    DynamicList<label> BladeP2ID(BladeP2.size());
                    for (int iID = 0; iID < BladeP2[m].size(); iID++)    
                    {
                        label a(iID);
                        BladeP2ID.append(a);
                    }
//                 Info << "BladeP2ID = " << BladeP2ID << endl;
                    List<label> BladeP2IDMP = 
                        extrapolate
                        (
                                    TC,
                                    BladeP2[m],
                                    BladeP2ID
                        );
//                 Info << "BladeP2IDMP = " << BladeP2IDMP << endl;
                    DynamicList<scalar> clTemp(2);
                    DynamicList<scalar> cdTemp(2);
                    DynamicList<scalar> BladeP2Temp(2);
                    for ( int iTemp = 0; iTemp < 2; iTemp++ )
                    {
                        clTemp[iTemp] = 
                                interpolate
                                    (
                                        alphaC[i][j][k],
                                        airfoilAlpha[airfoilsMP[iTemp]],
                                        airfoilCl[airfoilsMP[iTemp]]
                                    );
                        cdTemp[iTemp] = 
                                interpolate
                                    (
                                        alphaC[i][j][k],
                                        airfoilAlpha[airfoilsMP[iTemp]],
                                        airfoilCd[airfoilsMP[iTemp]]
                                    );
                        BladeP2Temp[iTemp] = BladeP2[m][BladeP2IDMP[iTemp]];
                    }
//                 Info << "TC = " << TC << endl;
//                 Info << "BladeP2Temp \t clTemp \t cdTemp " << endl;
//                 for ( int iTemp = 0; iTemp < 2; iTemp++ )
//                 {
//                 Info << BladeP2Temp[iTemp] << "\t" << clTemp[iTemp] << "\t" << cdTemp[iTemp] << endl;
//                 }
                    if (BladeP2Temp[0] == BladeP2Temp[1])
                    {
                        Cl[i][j][k] = clTemp[0];
                        Cd[i][j][k] = cdTemp[0];
                    }
                    else
                    {
                        Cl[i][j][k] = 
                        interpolate
                        (    
                            TC,
                            BladeP2Temp,
                            clTemp
                        );
                        Cd[i][j][k] = 
                        interpolate
                        (    
                            TC,
                            BladeP2Temp,
                            cdTemp
                        );
                    }
//                 if(i == 0 && j == 0)
//                 {
//                 Info << bladeRadius[i][j][k] << "\t" << TC << "\t" << Cl[i][j][k] << "\t" << Cd[i][j][k] << endl;
//                 Info << "Cl = " << Cl[i][j][k] << endl;
//                 Info << "Cd = " << Cd[i][j][k] << endl;
//                 }
                }
                else
                {
					FatalErrorInFunction
		        	        << " Invalid airfoil interpolation method "<<airfoilInterpolationMethod[m]<< "." 
		        	        << " Valid sampling methods are \n" 
		        	        << " LinearInterpolation \n" 
		        	        << " IDInterpolation \n" 
		        	        << " TuInterpolation \n"  
		        	        << " TCInterpolation \n"
		        	        << exit(FatalError);

                }


                // Calculate the local epsilon parameter and smearing radius
                //- Take global (constant) values along the blade
                if (smearingModel[i] == "Global")
                {
                    epsilonC[i][j][k] = epsilon[i];
                    smearingRadiusC[i][j][k] = smearRadius[i];
                }

                // - Local chord based value of epsilon
                else if (smearingModel[i] == "localChord")
                {
                    epsilonC[i][j][k] = localChordFraction[i]*Chord[i][j][k];
                    smearingRadiusC[i][j][k] = 3.0*epsilonC[i][j][k]; 
                }

                // - Lift Coefficient based model of epsilon  
                else if (smearingModel[i] == "liftCoefficient")
                {
                    epsilonC[i][j][k] = localChordFraction[i]*Chord[i][j][k]*Cl[i][j][k];
                    smearingRadiusC[i][j][k] = 3.0*epsilonC[i][j][k];                  
                }

                // - Local chord based value of epsilon plus a cap determined by input
                else if (smearingModel[i] == "localChordCapped")
                {
                    epsilonC[i][j][k] = min(localChordFraction[i]*Chord[i][j][k],epsilon[i]);
                    smearingRadiusC[i][j][k] = 3.0*epsilonC[i][j][k]; 
                }

                // - Local chord based value of epsilon plus a cap determined by input
                else if (smearingModel[i] == "powerFunction")
                {
                    epsilonC[i][j][k] =      epsilon[i] +
                                       -1.0*(epsilon[i] - localChordFraction[i]*Chord[i][j][Chord[i][j].size()-1]) * 
                                        (pow((bladeRadius[i][j][k]-HubRad[i])/(TipRad[m]-HubRad[i]), epsilonFunctionPow[i]));
                    smearingRadiusC[i][j][k] = 3.0*epsilonC[i][j][k]; 
                }

                // - Local chord based value of epsilon plus a cap determined by input
                else if (smearingModel[i] == "cosineFunction")
                {
                    epsilonC[i][j][k] = 0.5*(epsilon[i] + localChordFraction[i]*Chord[i][j][Chord[i][j].size()-1]) +
                                        0.5*(epsilon[i] - localChordFraction[i]*Chord[i][j][Chord[i][j].size()-1]) * 
                                        cos((Foam::constant::mathematical::pi)*pow((bladeRadius[i][j][k]-HubRad[i])/(TipRad[m]-HubRad[i]), epsilonFunctionPow[i]));
                    smearingRadiusC[i][j][k] = 3.0*epsilonC[i][j][k]; 
                }
                else
                {
					FatalErrorInFunction
		        	        << " Invalid Smearing Model "<<smearingModel[i]<< "." 
		        	        << " Valid smearing models are \n" 
		        	        << " Global \n" 
		        	        << " localChord \n" 
		        	        << " liftCoefficient \n"  
		        	        << " localChordCapped \n"
		        	        << " powerFunction \n"
		        	        << " cosineFunction \n"
		        	        << exit(FatalError);
                }

                // Apply tip/root-loss correction factor.
                scalar F    = 1.0;
                scalar Faxi = 1.0;
                scalar Ftan = 1.0;
                scalar Cl_add = 0.0;
                scalar Cd_add = 0.0;

                if (tipRootLossCorrType[i] == "none")
                {
                    F = 1.0;
                }

                else if (tipRootLossCorrType[i] == "None")
                {
                    F = 1.0;
                }

                else if (tipRootLossCorrType[i] == "Glauert")
                {
                    // Glauert correction adopted from NREL SOWFA code
                    scalar windAng = alphaC[i][j][k] + (Twist[i][j][k] + pitch[i]); // deg

                    scalar g = 1.0;

                    scalar ftip  = (TipRad[m] - bladeRadius[i][j][k]) / ((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar argTip1 = -g * (NumBl[m] / 2.0) * ftip;
                    scalar argTip2 = exp(max(min(30, argTip1),-30));
                    scalar Ftip = (2.0/(Foam::constant::mathematical::pi)) * acos(max(min(1.0,argTip2),-1.0));

                    scalar froot = (bladeRadius[i][j][k] - HubRad[i]) / ((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar argRoot1 = -g * (NumBl[m] / 2.0) * froot;
                    scalar argRoot2 = exp(max(min(30, argRoot1),-30));
                    scalar Froot = (2.0/(Foam::constant::mathematical::pi)) * acos(max(min(1.0,argRoot2),-1.0));

                    F = Ftip * Froot;
                }
                else if (tipRootLossCorrType[i] == "Shen")
                {
                    // Shen correction 'Shen et al. (2005) Tip Loss 
                    // Correction for Wind Turbine Computations, 
                    // Wind Energy 8:457-475'
                    scalar windAng = alphaC[i][j][k] + (Twist[i][j][k] + pitch[i]); // deg

                    scalar g = exp(-0.125*((NumBl[m] * tipSpeedRatioShen[i])-21.0)) + 0.1;

                    scalar ftip = (TipRad[m] - bladeRadius[i][j][k]) / ((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar arg1 = -g * (NumBl[m] / 2.0) * ftip; 
                    scalar arg2 = exp(max(min(30, arg1),-30));
                    scalar Ftip = (2.0/(Foam::constant::mathematical::pi))* 
                                        acos(max(min(1.0,arg2),-1.0));
                    F = Ftip;
                }
                else if (tipRootLossCorrType[i] == "Wimshurst")
                {
                    scalar windAng = alphaC[i][j][k] + (Twist[i][j][k] + pitch[i]); // deg

//                  scalar c1_ax = 0.1219; // Wimshurst & Willden 2017
//                  scalar c2_ax = 21.52;  // Wimshurst & Willden 2017
                    scalar g_axi = exp(-c1Faxi[i]*((NumBl[m] * tipSpeedRatioShen[i])-c2Faxi[i])) + c3Faxi[i];
                    scalar ftip_axi = (TipRad[m] - bladeRadius[i][j][k])/((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar arg1_axi = -g_axi * (NumBl[m] / 2.0) * ftip_axi;
                    scalar arg2_axi = exp(max(min(30, arg1_axi),-30));
                    scalar Ftip_axi = (2.0/(Foam::constant::mathematical::pi)) * acos(max(min(1.0,arg2_axi),-1.0));
                    Faxi = Ftip_axi;

        
//                  scalar c1_tan = 0.0984; // Wimshurst & Willden 2017
//                  scalar c2_tan = 13.026; // Wimshurst & Willden 2017
                    scalar g_tan = exp(-c1Ftan[i]*((NumBl[m] * tipSpeedRatioShen[i])-c2Ftan[i])) + c3Ftan[i];
                    scalar ftip_tan  = (TipRad[m] - bladeRadius[i][j][k]) / ((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar arg1_tan = -g_tan * (NumBl[m] / 2.0) * ftip_tan; 
                    scalar arg2_tan = exp(max(min(30, arg1_tan),-30));
                    scalar Ftip_tan = (2.0/(Foam::constant::mathematical::pi)) * acos(max(min(1.0,arg2_tan),-1.0));
                    Ftan = Ftip_tan;
                }
                else if (tipRootLossCorrType[i] == "CaoWillden") // Cao & Willden 2020
                {
                    // Glauert correction adopted from NREL SOWFA code
                    scalar windAng = alphaC[i][j][k] + (Twist[i][j][k] + pitch[i]); // deg

                    scalar g = 1.0;

                    scalar ftip  = (TipRad[m] - bladeRadius[i][j][k]) / ((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar argTip1 = -g * (NumBl[m] / 2.0) * ftip;
                    scalar argTip2 = exp(max(min(30, argTip1),-30));
                    scalar Ftip = (2.0/(Foam::constant::mathematical::pi)) * acos(max(min(1.0,argTip2),-1.0));

                    scalar froot = (bladeRadius[i][j][k] - HubRad[i]) / ((bladeRadius[i][j][k] * sin(windAng*degRad))+SMALL);
                    scalar argRoot1 = -g * (NumBl[m] / 2.0) * froot;
                    scalar argRoot2 = exp(max(min(30, argRoot1),-30));
                    scalar Froot = (2.0/(Foam::constant::mathematical::pi)) * acos(max(min(1.0,argRoot2),-1.0));

                    F = Ftip * Froot;

                    // Cao & Willden 2020 tip correction added to Glauert
                    scalar gL = -1.9;
                    scalar gD =  0.1;

                    p_gamma += bladeRadius[i][j][k] * rotSpeed[i] * (0.5*F*Cl[i][j][k]*VmagC[i][j][k]*VmagC[i][j][k]*Chord[i][j][k])/(2*Chord[i][j][k]*VmagC[i][j][k]);

                    Cl_add = gL * (p_gamma / (0.5*max(1.0,VmagC[i][j][k])*max(1.0,VmagC[i][j][k])));
                    Cd_add = gD * (p_gamma / (0.5*max(1.0,VmagC[i][j][k])*max(1.0,VmagC[i][j][k])));
                }
				else
				{
					FatalErrorInFunction
		                << " Invalid tip correction method "<<tipRootLossCorrType[i]<< "." 
		                << " Valid tip corrections are \n" 
		                << " None \n" 
		                << " Glauert \n"  
		                << " Shen \n"
		                << " Wimshurst \n"
		                << " CaoWillden \n"
		                << exit(FatalError);
				}

                //- Apply a dynamic gust model 
                scalar ClUnsteady = 0.0;

                //- Placeholder for when this is implemented
                if (dynamicGustModel[i] == "None")
                {
                    ClUnsteady =  0.0;
                }

                //- Placeholder for when this is implemented
                else if (dynamicGustModel[i] == "none")
                {
                    ClUnsteady =  0.0;
                }

                //- Placeholder for when this is implemented
                else if (dynamicGustModel[i] == "Theodorsen")
                {
                    ClUnsteady = 0.0;
                }

                //- Add the dynamic gust contribution to the quasi-static lift coefficient
                Cl[i][j][k] += ClUnsteady;

                // Check Cl and Cd for Tu Interpolation
//              if(i==0 and j==0) {Info << Cl[i][j][k] << "\t" << Cd[i][j][k] << endl;}

                // Using Cl, Cd, wind velocity, chord, and actuator element
                // width, calculate the lift and drag per density.
                lift[i][j][k] =
                    0.5*F*(Cl[i][j][k]+Cl_add)*VmagC[i][j][k]*VmagC[i][j][k]*Chord[i][j][k]*db[i][k];
                drag[i][j][k] =
                    0.5*F*(Cd[i][j][k]+Cd_add)*VmagC[i][j][k]*VmagC[i][j][k]*Chord[i][j][k]*db[i][k];

                // Rotate chordwise chordAlignedVector about (its) z axis by
                // effective AoA
                vector axis = chordAlignedVectors[i][j][k][2];
                // clockwise rotation about z-chordAlignedVector
                scalar angle = -alphaC[i][j][k]*degRad;
                vector dragVector =
                    rotatePoint
                    (
                        chordAlignedVectors[i][j][k][1],
                        vector::zero,
                        axis,
                        angle
                    );

                dragVector = dragVector/mag(dragVector);

                vector liftVector = dragVector^bladeAlignedVectors[i][j][2];
                liftVector = liftVector/mag(liftVector);

                // reversed to give blade-on-fluid force:
                liftVector = -lift[i][j][k] * liftVector;
                dragVector = -drag[i][j][k] * dragVector;

                // Add up lift and drag to get the resultant force/density
                // applied to this blade element.
                bladeForce[i][j][k] = liftVector + dragVector;
//              Info<<"Check: "<<bladeForce[i][j][k]<<endl;

                // Find the component of the blade element force/density in the
                // axial (along the shaft) direction.
                axialForce[i][j][k] = Faxi*(-bladeForce[i][j][k] & uvShaft[i]);

                // Find the component of the blade element force/density in the
                // tangential (torque-creating) direction.
                tangentialForce[i][j][k] =
                    Ftan*(bladeForce[i][j][k] & bladeAlignedVectors[i][j][1]);
        
        
                bladeForce[i][j][k] = -axialForce[i][j][k]*uvShaft[i] + tangentialForce[i][j][k]*bladeAlignedVectors[i][j][1];
//              Info<<"Check: "<<bladeForce[i][j][k]<<endl;
//              Info<<"Check: Axial Norm= "<<axialForce[i][j][k]<<",dir: "<<uvShaft[i]<<endl;
//              Info<<"Check: Tangential Norm= "<<tangentialForce[i][j][k]<<",dir: "<<bladeAlignedVectors[i][j][1]<<endl;

                // Find the component of the blade element force/density in the
                // normal (to the chord) direction.
                normalForce[i][j][k] = -bladeForce[i][j][k] & chordAlignedVectors[i][j][k][0];

                // Find the component of the blade element force/density in the
                // chordwise (parallel to the chord) direction.
                chordwiseForce[i][j][k] = bladeForce[i][j][k] & chordAlignedVectors[i][j][k][1];

                // Add this blade element's contribution to thrust to the total
                // turbine thrust.
                thrust[i] += axialForce[i][j][k];

                // Add this blade element's contribution to torque to the total
                // turbine torque.
                torqueRotor[i] +=
                    tangentialForce[i][j][k]*bladeRadius[i][j][k]
                   *cos(PreCone[m][j]);
            }
        }
        // Compute power based on torque and rotation speed.
        powerRotor[i] = torqueRotor[i] * rotSpeed[i];
    }
}

void Foam::fv::actuationLineSource::computeBladeForceCorrector()
{

    DynamicList<List<List<scalar> > > cellIntegration;

    // Proceed turbine by turbine.
    forAll(bladeForceCorrector, i)
    {
        // Initialise the correction function
        cellIntegration.append
        (    
            List<List<scalar> >(NumBl[i], List<scalar>(numBladePoints[i],0.0))
        );

        if (bladeForceCorrector_[i] == 1)
        {
            // Proceed over all sphere cells of that turbine if there are sphere
            // cells.
            if (allSphereCellsI[i][Pstream::myProcNo()].size() > 0)
            {
                forAll(allSphereCellsI[i][Pstream::myProcNo()], m)
                {
    
                    // For each blade.
                    forAll(bladeForceCorrector[i], j)
                    {
    
                        // For each blade point.
                        forAll(bladeForceCorrector[i][j], k)
                        {
                            label proc = Pstream::myProcNo();
                            scalar dis =
                                mag(mesh_.C()[allSphereCellsI[i][proc][m]]
                              - bladePoints[i][j][k]);
                            if (dis <= smearingRadiusC[i][j][k])
                            {
                                cellFinder[allSphereCellsI[i][proc][m]] += 1.0;
                                cellIntegration[i][j][k] += 
                                    mesh_.V()[allSphereCellsI[i][proc][m]]
                                   *(
                                        Foam::exp(-Foam::sqr(dis/epsilonC[i][j][k]))
                                       /(
                                            Foam::pow(epsilonC[i][j][k],3)*Foam::pow
                                            (
                                                Foam::constant::mathematical::pi,1.5
                                            )
                                        )
                                    );
                            }
                        }
                    }
                }
            }

            forAll(cellIntegration[i],j)
            {
                forAll(cellIntegration[i][j],k)
                {
                    reduce(cellIntegration[i][j][k],sumOp<scalar>());
                }
            }

            forAll(cellIntegration[i],j)
            {
                forAll(cellIntegration[i][j],k)
                {
                    bladeForceCorrector[i][j][k] = 1.0 / max(cellIntegration[i][j][k],1.0E-2);
                }
            }
        }
        else
        {
            forAll(bladeForceCorrector[i],j)
            {
                forAll(bladeForceCorrector[i][j],k)
                {
                    bladeForceCorrector[i][j][k] = 1.0;
                }
            }
        }
    }
}


void Foam::fv::actuationLineSource::computeBodyForce()
{
    // Zero out the body force to begin with.
    bodyForce *= 0.0;

    // Proceed turbine by turbine.
    forAll(bladeForce, i)
    {

        // Proceed over all sphere cells of that turbine if there are sphere
        // cells.
        if (allSphereCellsI[i][Pstream::myProcNo()].size() > 0)
        {
            forAll(allSphereCellsI[i][Pstream::myProcNo()], m)
            {

                // For each blade.
                forAll(bladeForce[i], j)
                {

                    // For each blade point.
                    forAll(bladeForce[i][j], k)
                    {
                        label proc = Pstream::myProcNo();
                        scalar dis =
                            mag(mesh_.C()[allSphereCellsI[i][proc][m]]
                          - bladePoints[i][j][k]);
                        if (dis <= smearingRadiusC[i][j][k])
                        {
                            cellFinder[allSphereCellsI[i][proc][m]] += 1.0;
                            bodyForce[allSphereCellsI[i][proc][m]] +=
                                bladeForce[i][j][k]*bladeForceCorrector[i][j][k]
                               *(
                                    Foam::exp(-Foam::sqr(dis/epsilonC[i][j][k]))
                                   /(
                                        Foam::pow(epsilonC[i][j][k],3)*Foam::pow
                                        (
                                            Foam::constant::mathematical::pi,1.5
                                        )
                                    )
                                );
                        }
                    }
                }
            }
        }
    }
}

vector
Foam::fv::actuationLineSource::rotatePoint
(
    vector point,
    vector rotationPoint,
    vector ax,
    scalar angle
)
{
    // Declare and define the rotation matrix.
    tensor RM;
    RM.xx() = Foam::sqr(ax.x()) + (1.0 - Foam::sqr(ax.x())) * Foam::cos(angle);
    RM.xy() = ax.x()*ax.y()*(1.0 - Foam::cos(angle)) - ax.z()*Foam::sin(angle);
    RM.xz() = ax.x()*ax.z()*(1.0 - Foam::cos(angle)) + ax.y()*Foam::sin(angle);
    RM.yx() = ax.x()*ax.y()*(1.0 - Foam::cos(angle)) + ax.z()*Foam::sin(angle);
    RM.yy() = Foam::sqr(ax.y()) + (1.0 - Foam::sqr(ax.y()))*Foam::cos(angle);
    RM.yz() = ax.y()*ax.z()*(1.0 - Foam::cos(angle)) - ax.x()*Foam::sin(angle);
    RM.zx() = ax.x()*ax.z()*(1.0 - Foam::cos(angle)) - ax.y()*Foam::sin(angle);
    RM.zy() = ax.y()*ax.z()*(1.0 - Foam::cos(angle)) + ax.x()*Foam::sin(angle);
    RM.zz() = Foam::sqr(ax.z()) + (1.0 - Foam::sqr(ax.z()))*Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    point = point - rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point.
    point = point + rotationPoint;

    return point;
}


scalar
Foam::fv::actuationLineSource::interpolate
(
    scalar xNew,
    DynamicList<scalar>& xOld,
    DynamicList<scalar>& yOld
)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return
            yOld[indexM]
          + (
                (yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM])
            )
           *(xNew - xOld[indexM]);
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return
            yOld[indexM]
          + (
                (yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM])
            )
           *(xNew - xOld[indexM]);
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}

List<label>
Foam::fv::actuationLineSource::extrapolate
(
    scalar xNew,
    DynamicList<scalar>& xOld,
    DynamicList<label>& yOld
)
{
    label indexP = 0;
    label indexM = 0;
    List<label> labelsMP(2,0);
    if( xOld[0] < xOld[1] )
    {
        for (int i = 0; i < xOld.size()-1; i++ )
        {
            if ( (xNew > xOld[i]) && (xNew < xOld[i+1] ) )
            {
            indexM = i;
            break;
            }
        }
        indexP = indexM + 1;
        if ( xNew > xOld[xOld.size()-1] ) 
        {
            indexP = yOld[xOld.size()-1]; 
            indexM = indexP;
        }
        if ( xNew < xOld[0] ) 
        {
            indexP = yOld[0]; 
            indexM = indexP;
        }
    }
    else if ( xOld[0] > xOld[1] )
    {
        for (int i = 0; i < xOld.size()-1; i++ )
        {
            if ( (xNew < xOld[i]) && (xNew > xOld[i+1] ) )
            {
            indexM = i;
            break;
            }
        }
        indexP = indexM + 1;
        if ( xNew < xOld[xOld.size()-1] ) 
        {
            indexP = yOld[xOld.size()-1]; 
            indexM = indexP;
        }
        if ( xNew > xOld[0] ) 
        {
            indexP = yOld[0]; 
            indexM = indexP;
        }
    }

    labelsMP[0] = yOld[indexM];
    labelsMP[1] = yOld[indexP];
    return labelsMP;
}

label
Foam::fv::actuationLineSource::interpolate
(
    scalar xNew,
    DynamicList<scalar>& xOld,
    DynamicList<label>& yOld
)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return
            round
            (
                yOld[indexM]
              + (
                    (yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM])
                )
               *(xNew - xOld[indexM])
            );
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return
            round
            (
                yOld[indexM]
              + (
                    (yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM])
                )
               *(xNew - xOld[indexM])
            );
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}


void Foam::fv::actuationLineSource::openOutputFiles()
{
    if (Pstream::master())
    {
        // Create the name of the root of where turbine files get ouput.
        fileName rootDir;

        if (Pstream::parRun())
        {
            rootDir = runTime_.path()/"../turbineOutput";
        }
        else
        {
            rootDir = runTime_.path()/"turbineOutput";
        }

        // Check to see if the turbineOutput directory exists; if not, create
        // it.
        if (!isDir(rootDir))
        {
            mkDir(rootDir);
        }

        // Check to see if the start time directory exists within the
        // turbineOutput directory; if not, create it.
        if (!isDir(rootDir/time))
        {
            mkDir(rootDir/time);
        }

        // Write a file with turbineArrayProperties setups
        setupFile_ = new OFstream(rootDir/time/"turbineArrayProperties");
        forAll(bladePoints,i)
        {
            *setupFile_
                << "turbine" << i << endl;
            *setupFile_
                << "{" << endl;
            *setupFile_
                << "    turbineType" << "    " << turbineType[i] << endl;
            *setupFile_
                << "    baseLocation" << "    " << baseLocation[i] << endl;
            *setupFile_
                << "    numBladePoints" << "    " << numBladePoints[i] << endl;
            *setupFile_
                << "    pointDistType" << "    " << pointDistType[i] << endl;
            *setupFile_
                << "    bladeForceCorrector" << "    " << bladeForceCorrector_[i] << endl;
            *setupFile_
                << "    samplingMethod" << "    " << samplingMethod[i] << endl;
            *setupFile_
                << "        sampleDistScalar" << "    " << sampleDistScalar[i] << endl;
            *setupFile_
                << "        sampleDistParam" << "    " << sampleDistParam[i] << endl;
            *setupFile_
                << "        numSamplePoints" << "    " << numSamplePoints[i] << endl;
            *setupFile_
                << "    smearingModel" << "    " << smearingModel[i] << endl;
            *setupFile_
                << "        epsilon" << "    " << epsilon[i] << endl;
            *setupFile_
                << "        localChordFraction" << "    " << localChordFraction[i] << endl;
            *setupFile_
                << "        epsilonFunctionPow" << "    " << epsilonFunctionPow[i] << endl;
            *setupFile_
                << "    sphereRadiusScalar" << "    " << sphereRadiusScalar[i] << endl;
            *setupFile_
                << "    tipRootLossCorrType" << "    " << tipRootLossCorrType[i] << endl;
            *setupFile_
                << "    tipSpeedRatioShen" << "    " << tipSpeedRatioShen[i] << endl;
            *setupFile_
                << "        c1Faxi" << "    " << c1Faxi[i] << endl;
            *setupFile_
                << "        c2Faxi" << "    " << c2Faxi[i] << endl;
            *setupFile_
                << "        c3Faxi" << "    " << c3Faxi[i] << endl;
            *setupFile_
                << "        c1Ftan" << "    " << c1Ftan[i] << endl;
            *setupFile_
                << "        c2Ftan" << "    " << c2Ftan[i] << endl;
            *setupFile_
                << "        c3Ftan" << "    " << c3Ftan[i] << endl;
            *setupFile_
                << "    airfoilInterpolationMethod" << "    " << airfoilInterpolationMethod[i] << endl;
            *setupFile_
                << "    dynamicGustModel" << "    " << dynamicGustModel[i] << endl;
            *setupFile_
                << "    rotationDir" << "    " << rotationDir[i] << endl;
            *setupFile_
                << "    Azimuth" << "    " << azimuth[i] << endl;
            *setupFile_
                << "    RotSpeed" << "    " << rotSpeed[i] << endl;
            *setupFile_
                << "    NacYaw" << "    " << nacYaw[i] << endl;
            *setupFile_
                << "    fluidDensity" << "    " << fluidDensity[i] << endl;
            *setupFile_
                << "    rotorUref" << "    " << rotorUref[i] << endl;
            *setupFile_
                << "    rotorUturb" << "    " << rotorUturb[i] << endl;
            *setupFile_
                << "    pitchZero" << "    " << pitchZero[i] << endl;
            *setupFile_
                << "}" << endl;
        }

        // Create a total aerodynamic torque file.
        //torqueRotorFile(rootDir/time/"torqueRotor");
        torqueRotorFile_ = new OFstream(rootDir/time/"torqueRotor");
        *torqueRotorFile_
            << "#Turbine    Time(s)    dt(s)    rotor torque (N-m)" << endl;

        // Create a total thrust file.
        thrustFile_ = new OFstream(rootDir/time/"thrust");
        *thrustFile_
            << "#Turbine    Time(s)    dt(s)    thrust (N)" << endl;

        // Create a total power file.
        powerRotorFile_ = new OFstream(rootDir/time/"powerRotor");
        *powerRotorFile_
            << "#Turbine    Time(s)    dt(s)    rotor power (W)" << endl;


        // Create a blade points file.
        radiusCFile_ = new OFstream(rootDir/time/"radiusC");
        *radiusCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    radiusC(m)" << endl;
        forAll (bladePoints,i)
        {
            forAll(bladePoints[i],j)
            {
                *radiusCFile_<< i << " " << j << " " ; 
                forAll(bladePoints[i][j],k)
                {
                    *radiusCFile_<<mag(bladePoints[i][j][k]-rotorApex[i])<<" ";
                }
                *radiusCFile_ << endl;
            }
            *radiusCFile_ << endl;
        }


        // Create a rs file.
        rsFile_ = new OFstream(rootDir/time/"rs");
        *rsFile_
            << "#Turbine    Blade    Time(s)    dt(s)    rs(m)" << endl;

        // Create a qC file.
        qCFile_ = new OFstream(rootDir/time/"qC");
        *qCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    qC(m/s)" << endl;

        // Create a Guassian smearing parameter file.
        epsilonCFile_ = new OFstream(rootDir/time/"epsilonC");
        *epsilonCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    epsilon(m)" << endl;

        // Create a smearing radius file.
        smearingRadiusCFile_ = new OFstream(rootDir/time/"smearingRadiusC");
        *smearingRadiusCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    smearing radius R_s(m)" << endl;

        // Create a circulation strength file.
        gammaCFile_ = new OFstream(rootDir/time/"gammaC");
        *gammaCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    Circulation Strength Gamma_C(m2/s)" << endl;

        // Create effective velocity magnitude file.
        VmagCFile_ = new OFstream(rootDir/time/"VmagC");
        *VmagCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    VmagC(m/s)" << endl;

        // Create an angle of attack file.
        alphaCFile_ = new OFstream(rootDir/time/"alphaC");
        *alphaCFile_
            << "#Turbine    Blade    Time(s)    dt(s)    angle-of-attack(degrees)"
            << endl;

        // Create a coefficient of lift file.
        ClFile_ = new OFstream(rootDir/time/"Cl");
        *ClFile_
            << "#Turbine    Blade    Time(s)    dt(s)    Cl" << endl;

        // Create a coefficient of drag file.
        CdFile_ = new OFstream(rootDir/time/"Cd");
        *CdFile_
            << "#Turbine    Blade    Time(s)    dt(s)    Cd" << endl;

        // Create a lift file.
        liftFile_ = new OFstream(rootDir/time/"lift");
        *liftFile_
            << "#Turbine    Blade    Time(s)    dt(s)    lift (N)" << endl;

        // Create a drag file.
        dragFile_ = new OFstream(rootDir/time/"drag");
        *dragFile_
            << "#Turbine    Blade    Time(s)    dt(s)    drag (N)" << endl;

        // Create a axial force file.
        axialForceFile_ = new OFstream(rootDir/time/"axialForce");
        *axialForceFile_
            << "#Turbine    Blade    Time(s)    dt(s)    axial force (N)"
            << endl;

        // Create a tangential force file.
        tangentialForceFile_ = new OFstream(rootDir/time/"tangentialForce");
        *tangentialForceFile_
            << "#Turbine    Blade    Time(s)    dt(s)    tangential force (N)"
            << endl;

    // mark additions {
        // Create a normal force file.
        normalForceFile_ = new OFstream(rootDir/time/"normalForce");
        *normalForceFile_
            << "#Turbine    Blade    Time(s)    dt(s)    normal force (N)"
            << endl;

        // Create a chordwise force file.
        chordwiseForceFile_ = new OFstream(rootDir/time/"chordwiseForce");
        *chordwiseForceFile_
            << "#Turbine    Blade    Time(s)    dt(s)    chordwise force (N)"
            << endl;

        // Create a bladeForceCorrector check file.
        bladeForceCorrectorFile_ = new OFstream(rootDir/time/"bladeForceCorrector");
        *bladeForceCorrectorFile_
            << "#Turbine    Blade    Time(s)    dt(s)    locations (m)" << endl;


        // Create a samplePoints location check file.
//      samplePointsFile_ = new OFstream(rootDir/time/"samplePointsFile");
//      *samplePointsFile_
//          << "#Turbine    Blade    Time(s)    dt(s)    locations (m)" << endl;

//      // Create a probe point check file.
//      probeCellIndexFile_ = new OFstream(rootDir/time/"probeCellIndexFile");
//      *probeCellIndexFile_
//          << "#Turbine    Blade    Time(s)    dt(s)    pointIndex" << endl;

//      // Create a samplePoints check file.
//      samplePointsUFile_ = new OFstream(rootDir/time/"samplePointsUFile");
//      *samplePointsUFile_
//          << "#Turbine    Blade    Time(s)    dt(s)    U (m/s)" << endl;
     //mark additions }
    }
}


void Foam::fv::actuationLineSource::printOutputFiles()
{
    if (Pstream::master())
    {
        forAll(bladePoints,i)
        {
            // Write out time and delta t.
            *torqueRotorFile_ << i << " " << time << " " << dt << " ";
            *thrustFile_ << i << " " << time << " " << dt << " ";
            *powerRotorFile_ << i << " " << time << " " << dt << " ";

            // Write out information for each turbine.
            *torqueRotorFile_ << torqueRotor[i]*fluidDensity[i] << endl;
            *thrustFile_ << thrust[i]*fluidDensity[i] << endl;
            *powerRotorFile_ << powerRotor[i]*fluidDensity[i] << endl;

            // Proceed blade by blade.
            forAll(bladePoints[i], j)
            {
                // Write out time and delta t.
                
                *rsFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *qCFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *VmagCFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *alphaCFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *epsilonCFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *smearingRadiusCFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *gammaCFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *ClFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *CdFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *liftFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *dragFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *axialForceFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *tangentialForceFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *normalForceFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *chordwiseForceFile_
                    << i << " " << j << " " << time << " " << dt << " ";
                *bladeForceCorrectorFile_
                    << i << " " << j << " " << time << " " << dt << " ";

                if(i == 0 && j== 0)
                {
//                  *bladePointsFile_
//                      << i << " " << j << " " << time << " " << dt;
//                  *samplePointsFile_
//                      << i << " " << j << " " << time << " " << dt;
//                  *probeCellIndexFile_
//                      << i << " " << j << " " << time << " " << dt;
//                  *samplePointsUFile_
//                      << i << " " << j << " " << time << " " << dt;
                }

                forAll(bladePoints[i][j], k)
                {

                    *rsFile_
                        << rs[i][j][k] << " ";
                    *qCFile_
                        << qC[i][j][k] << " ";
                    *VmagCFile_
                        << VmagC[i][j][k] << " ";
                    *alphaCFile_
                        << alphaC[i][j][k] << " ";
                    *epsilonCFile_
                        << epsilonC[i][j][k] << " ";
                    *smearingRadiusCFile_
                        << smearingRadiusC[i][j][k] << " ";
                    *gammaCFile_
                        << GammaC[i][j][k] << " ";
                    *ClFile_
                        << Cl[i][j][k] << " ";
                    *CdFile_
                        << Cd[i][j][k] << " ";
                    *liftFile_
                        << lift[i][j][k]*fluidDensity[i] << " ";
                    *dragFile_
                        << drag[i][j][k]*fluidDensity[i] << " ";
                    *axialForceFile_
                        << axialForce[i][j][k]*fluidDensity[i] << " ";
                    *tangentialForceFile_
                        << tangentialForce[i][j][k]*fluidDensity[i] << " ";
                    *normalForceFile_
                        << normalForce[i][j][k]*fluidDensity[i] << " ";
                    *chordwiseForceFile_
                        << chordwiseForce[i][j][k]*fluidDensity[i] << " ";
                    *bladeForceCorrectorFile_
                        << bladeForceCorrector[i][j][k] << " ";

//                  if(i == 0 && j== 0 && k == 10)
//                  {
//                      forAll(samplePoints[i][j][k], l)
//                      {
//                          *samplePointsFile_
//                              << " " << samplePoints[i][j][k][l];
//                      }
//                      forAll(samplePointsI[i][j][k], l)
//                      {
//                          *probeCellIndexFile_
//                              << " " << samplePointsI[i][j][k][l];
//                      }
//                      forAll(samplePointsU[i][j][k], l)
//                      {
//                          *samplePointsUFile_
//                              << " " << samplePointsU[i][j][k][l];
//                      }
//                      *samplePointsFile_ << endl;
//                      *probeCellIndexFile_ << endl;
//                      *samplePointsUFile_ << endl;
//                  }
                }
                *rsFile_ << endl;
                *qCFile_ << endl;
                *VmagCFile_ << endl;
                *alphaCFile_ << endl;
                *epsilonCFile_ << endl;
                *smearingRadiusCFile_ << endl;
                *gammaCFile_ << endl;
                *ClFile_ << endl;
                *CdFile_ << endl;
                *liftFile_ << endl;
                *dragFile_ << endl;
                *axialForceFile_ << endl;
                *tangentialForceFile_ << endl;
                *normalForceFile_ << endl;
                *chordwiseForceFile_ << endl;
                *bladeForceCorrectorFile_ << endl;
            }
        }

        //*torqueRotorFile_ << endl;
        //*thrustFile_ << endl;
        //*powerRotorFile_ << endl;

        *rsFile_ << endl;
        *qCFile_ << endl;
        *VmagCFile_ << endl;
        *alphaCFile_ << endl;
        *gammaCFile_ << endl;
        *epsilonCFile_ << endl;
        *smearingRadiusCFile_ << endl;
        *ClFile_ << endl;
        *CdFile_ << endl;
        *liftFile_ << endl;
        *dragFile_ << endl;
        *axialForceFile_ << endl;
        *tangentialForceFile_ << endl;
        *normalForceFile_ << endl;
        *chordwiseForceFile_ << endl;
        *bladeForceCorrectorFile_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuationLineSource::actuationLineSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),

    //- Degrees to radians conversion factor.
    degRad((Foam::constant::mathematical::pi)/180.0),

    //- Revolutions per minute to radians per second conversion factor.
    rpmRadSec(2.0*(Foam::constant::mathematical::pi)/60.0),

    //- Set the pointer to runTime
    runTime_(mesh_.time()),

    // Set the time step size.
    dt(runTime_.deltaT().value()),

    // Set the current simulation time.
    time(runTime_.timeName()),

    //- Velocity sampling.
    ms(mesh),
    interpolationDict(coeffs_.subDict("interpolationScheme")),

    //- debug.
    debug_(coeffs_.lookupOrDefault<bool>("debug",0)),

    // Initialize the body force.
    bodyForce
    (
        IOobject
        (
            "bodyForce",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimForce/dimVolume/dimDensity,vector::zero)
    ),

    // Initialize the cellFinder.
    cellFinder
    (
        IOobject
        (
            "cellFinder",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cellFinder",dimless,0.0)
    )
{

    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation line zone: "
        << this->name() << endl;

    // Define dictionary that defines the turbine array.
    IOdictionary turbineArrayProperties
    (
        IOobject
        (
            "turbineArrayProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read in the turbine array properties dictionary.  This is the uppermost
    // level dictionary that describes where the turbines are, what kind they
    // are, their initial state, and information about how the actuator line
    // method is applied to each turbine.
    {
        List<word> listTemp = turbineArrayProperties.toc();
        for (int i = 0; i < listTemp.size(); i++)
        {
            if (listTemp[i] != "globalProperties")
            {
                turbineName.append(listTemp[i]);
            }
        }
    }

    numTurbines = turbineName.size();

    outputControl = turbineArrayProperties.subDict("globalProperties")
       .lookupOrDefault<word>("outputControl","timeStep");
    outputInterval = turbineArrayProperties.subDict("globalProperties")
       .lookupOrDefault<scalar>("outputInterval",1);
    lastOutputTime = runTime_.startTime().value();
    outputIndex = 0;

    forAll(turbineName,i)
    {
        samplingMethod.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("samplingMethod")));
        numSamplePoints.append(int(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("numSamplePoints"))));
        sampleDistScalar.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("sampleDistScalar"))));
        sampleDistParam.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("sampleDistParam"))));
        bladeForceCorrector_.append(int(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("bladeForceCorrector"))));
        turbineType.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("turbineType")));
        baseLocation.append(vector(turbineArrayProperties
           .subDict(turbineName[i]).lookup("baseLocation")));
        numBladePoints.append(int(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("numBladePoints"))));
        pointDistType.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("pointDistType")));
        epsilon.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("epsilon"))));
        smearRadius.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("smearRadius"))));
        sphereRadiusScalar.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("sphereRadiusScalar"))));
        smearingModel.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("smearingModel")));
        localChordFraction.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("localChordFraction"))));           
        epsilonFunctionPow.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("epsilonFunctionPow"))));           
        tipSpeedRatioShen.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("tipSpeedRatioShen"))));
        tipRootLossCorrType.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("tipRootLossCorrType")));
        c1Faxi.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("c1Faxi"))));
        c2Faxi.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("c2Faxi"))));
        c3Faxi.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("c3Faxi"))));
        c1Ftan.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("c1Ftan"))));
        c2Ftan.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("c2Ftan"))));
        c3Ftan.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("c3Ftan"))));
        dynamicGustModel.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("dynamicGustModel")));
        airfoilInterpolationMethod.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("airfoilInterpolationMethod")));
        rotationDir.append(word(turbineArrayProperties
           .subDict(turbineName[i]).lookup("rotationDir")));
        rotSpeed.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("RotSpeed"))));
        azimuth.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("Azimuth"))));
        nacYaw.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("NacYaw"))));
        fluidDensity.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("fluidDensity"))));
        rotorUref.append(vector(turbineArrayProperties
           .subDict(turbineName[i]).lookup("rotorUref")));
        rotorUturb.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("rotorUturb"))));
        pitchZero.append(scalar(readScalar(turbineArrayProperties
           .subDict(turbineName[i]).lookup("pitchZero"))));
    lastOutputTime = runTime_.startTime().value();
    outputIndex = 0;
    }


    // Catalog the various types of turbines.  For example if three turbines are
    // GE 1.5 and two turbines are Siemens 2.3 machines, then assign the GE 1.5
    // an ID of 0 and the Siemens 2.3 an ID of 1.
    numTurbinesDistinct = 1;
    {
        turbineTypeDistinct.append(turbineType[0]);
        forAll(turbineType,i)
        {
            bool flag = false;
            for (int j = 0; j < numTurbinesDistinct; j++)
            {
                if (turbineType[i] == turbineTypeDistinct[j])
                {
                   flag = true;
                }
            }
            if (flag == false)
            {
                numTurbinesDistinct++;
                turbineTypeDistinct.append(turbineType[i]);
            }
        }
    }
    forAll(turbineType,i)
    {
        for (int j = 0; j < numTurbinesDistinct; j++)
        {
            if (turbineType[i] == turbineTypeDistinct[j])
            {
                turbineTypeID.append(j);
            }
        }
    }


    // For each distinct turbine, read in properties of that turbine from
    // separate dictionaries.

    for (int i = 0; i < numTurbinesDistinct; i++)
    {
        // Declare the turbineProperties dictionary for the ith turbine.
        IOdictionary turbineProperties
        (
            IOobject
            (
                turbineTypeDistinct[i],
                runTime_.constant(),
                "turbineProperties",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read in the data.
        NumBl.append(scalar(readScalar(turbineProperties.lookup("NumBl"))));
        TipRad.append(scalar(readScalar(turbineProperties.lookup("TipRad"))));
        HubRad.append(scalar(readScalar(turbineProperties.lookup("HubRad"))));
        UndSling.append(scalar(readScalar(turbineProperties.lookup("UndSling"))));
        OverHang.append(scalar(readScalar(turbineProperties.lookup("OverHang"))));
        TowerHt.append(scalar(readScalar(turbineProperties.lookup("TowerHt"))));
        Twr2Shft.append(scalar(readScalar(turbineProperties.lookup("Twr2Shft"))));
        ShftTilt.append(scalar(readScalar(turbineProperties.lookup("ShftTilt"))));
        PitchRate.append(scalar(readScalar(turbineProperties.lookup("PitchRate"))));
        PitchAmplitude.append(scalar(readScalar(turbineProperties.lookup("PitchAmplitude"))));
        PitchControllerType.append(word(turbineProperties.lookup("PitchControllerType")));

        List<scalar> preConeList = turbineProperties.lookup("PreCone");
        PreCone.append(preConeList);

        List<word> airfoilTypeList = turbineProperties.lookup("Airfoils");
        AirfoilType.append(airfoilTypeList);

        List<List<scalar> > bladeDataList = turbineProperties.lookup("BladeData");
        BladeData.append(bladeDataList);

        List<List<scalar> > rotorDesignList = turbineProperties.lookup("rotorDesign");
        rotorDesign.append(rotorDesignList);

        DynamicList<scalar> station;
        DynamicList<scalar> p1;
        DynamicList<scalar> p2;
        DynamicList<scalar> p3;
        DynamicList<scalar> p4;
        DynamicList<scalar> TC;
        DynamicList<label> id;

        forAll(BladeData[i], j)
        {
            station.append(BladeData[i][j][0]);
                 p1.append(BladeData[i][j][1]);
                 p2.append(BladeData[i][j][2]);
                 id.append(BladeData[i][j][3]);
        }

        BladeStation.append(station);
        BladeP1.append(p1);
        BladeP2.append(p2);
        BladeAirfoilTypeID.append(id);

        station.clear();
        p1.clear();
        p2.clear();
        id.clear();
        
        forAll(rotorDesign[i],j)
        {
            station.append(rotorDesign[i][j][0]);
                 p3.append(rotorDesign[i][j][1]);
                 p4.append(rotorDesign[i][j][2]);
                 TC.append(rotorDesign[i][j][3]);
        }
    
        DesignStation.append(station);
        DesignChord.append(p3);
        DesignTwist.append(p4);
        DesignTC.append(TC);
    
        station.clear();
        p3.clear();
        p4.clear();
        TC.clear();
    }


    // Catalog the various distinct types of airfoils used in the various
    // distinct types of turbines.
    int numAirfoilsDistinct = 1;
    {
        airfoilTypesDistinct.append(AirfoilType[0][0]);
        forAll(AirfoilType,i)
        {
            forAll(AirfoilType[i],j)
            {
                bool flag = false;
                for (int k = 0; k < numAirfoilsDistinct; k++)
                {
                    if (AirfoilType[i][j] == airfoilTypesDistinct[k])
                    {
                        flag = true;
                    }
                }
                if (flag == false)
                {
                    numAirfoilsDistinct++;
                    airfoilTypesDistinct.append(AirfoilType[i][j]);
                }
            }
        }
    }


    // Reassign airfoil type IDs to blades of each turbine based on the global
    // distinct list of airfoils.
    forAll(BladeAirfoilTypeID,i)
    {
        forAll(BladeAirfoilTypeID[i],j)
        {
            for (int k = 0; k < numAirfoilsDistinct; k++)
            {
                if
                (
                    AirfoilType[i][BladeAirfoilTypeID[i][j]] ==
                        airfoilTypesDistinct[k]
                )
                {
                    BladeAirfoilTypeID[i][j] = k;
                }
            }
        }
    }


    // For each distinct airfoil, read in the lift and drag versus angle
    // of attach data.
    for (int i = 0; i < numAirfoilsDistinct; i++)
    {
        // Declare the airfoilsProperties dictionary for the ith airfoil.
        IOdictionary airfoilProperties
        (
            IOobject
            (
                airfoilTypesDistinct[i],
                runTime_.constant(),
                "airfoilProperties",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read in the data.
        List<List<scalar> > airfoilDataList = airfoilProperties.lookup("airfoilData");
        airfoilData.append(airfoilDataList);

        DynamicList<scalar> alphaInt;
        DynamicList<scalar> ClInt;
        DynamicList<scalar> CdInt;

        forAll(airfoilData[i],j)
        {
            alphaInt.append(airfoilData[i][j][0]);
            ClInt.append(airfoilData[i][j][1]);
            CdInt.append(airfoilData[i][j][2]);
        }

        airfoilAlpha.append(alphaInt);
        airfoilCl.append(ClInt);
        airfoilCd.append(CdInt);

        alphaInt.clear();
        ClInt.clear();
        CdInt.clear();
    }


    // Convert nacelle yaw from cardinal directions (like on a compass)
    // to the normal convention of 0 degrees on the + x axis with
    // positive degrees in the counter-clockwise direction.
    forAll(nacYaw, i)
    {
        if (nacYaw[i] > 180.0)
        {
            nacYaw[i] = nacYaw[i] - 180.0;
        }
        else
        {
            nacYaw[i] = nacYaw[i] + 180.0;
        }
        nacYaw[i] = 90.0 - nacYaw[i];
        if (nacYaw[i] < 0.0)
        {
            nacYaw[i] = nacYaw[i] + 360.0;
        }
    }


    // Restart flag initialise
    if(runTime_.value() > runTime_.deltaT().value())
    {
        for (int i = 0; i< numTurbines; i++)
        {
            restart_flag.append(int(1));
        }
    }
    else
    {
        for (int i = 0; i< numTurbines; i++)
        {
            restart_flag.append(int(0));
        }
    }

    // Convert quantities in degrees into radians (dynamic lists
    // have to be done in loops).
    azimuth  = degRad * azimuth;
    rotSpeed = rpmRadSec * rotSpeed;
    nacYaw   = degRad * nacYaw;
    ShftTilt = degRad * ShftTilt;
    forAll(PreCone,i)
    {
        PreCone[i] = degRad * PreCone[i];
    }


    // Calculate tower shaft intersection and rotor apex locations. (The
    // i-index is at the turbine array level for each turbine and the j-
    // index is for each type of turbine--if all turbines are the same, j-
    // is always 0.)  The rotor apex is not yet rotated for initial yaw;
    // that is done below.
    for (int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];
        towerShaftIntersect.append(baseLocation[i]);
        towerShaftIntersect[i].z() =
            towerShaftIntersect[i].z() + TowerHt[j] + Twr2Shft[j];
        rotorApex.append(towerShaftIntersect[i]);
        rotorApex[i].x() =
            rotorApex[i].x()
          + ((OverHang[j] + UndSling[j])*Foam::cos(ShftTilt[j]));
        rotorApex[i].z() =
            rotorApex[i].z()
         +  (OverHang[j] + UndSling[j])*Foam::sin(ShftTilt[j]);
    }


    // Generating list of cells in sphere about each turbine i.
    // Labels are then gathered to the primary node.
    Info<< "find allSphereCells ... " << endl;
    for (int i = 0; i < numTurbines; i++) // each turbine i
    {
        List<DynamicList<label> > sphereCellsILocal
        (
            Pstream::nProcs(), DynamicList<label>()
        );
        scalar sphereRadius = 0.0;
        int j = turbineTypeID[i]; // distinct turbine type index j
        // one precone angle per blade (for this distinct turbine type).
        forAll(PreCone[j],k)
        {
            scalar sphereRadiusI =
                sphereRadiusScalar[i]*Foam::sqrt
                (
                    Foam::sqr((OverHang[j] + UndSling[j])
                  + TipRad[j]*Foam::sin(PreCone[j][k]))
                  + Foam::sqr(TipRad[j]*Foam::cos(PreCone[j][k]))
                );
            if (sphereRadiusI > sphereRadius)
            {
                sphereRadius = sphereRadiusI;
            }
        }
        Info<< "sphereRadius, turbine[" << i << "] = " << sphereRadius << endl;
        forAll(mesh_.cells(),cellI)
        {
            if (mag(mesh_.C()[cellI] - towerShaftIntersect[i]) <= sphereRadius)
            {
                //cellFinder[cellI] = 1.0;
                sphereCellsILocal[Pstream::myProcNo()].append(cellI);
            }
        }
        // [n][labels] n processors in charge of [labels] for turbine i.
        Pstream::gatherList(sphereCellsILocal);
        // each proc needs full allSphereCellsI for velocity sampling
        Pstream::scatterList(sphereCellsILocal);
        // appended tor each turbine leading to [i][n][labels].
        allSphereCellsI.append(sphereCellsILocal);
        sphereCellsILocal.clear();
    }
    Info<< "done" << endl;


    // Create the actuator line points (not yet rotated for initial nacelle
    // yaw or initial rotor azimuth. i-index is at array level, j-index is
    // for the type of turbine, k-index is for each blade, and m-index is
    // for each actuator point.  Also create other important vectors, and
    // initialize the blade force, blade aligned coordinate system, and
    // wind vectors to zero.
    totBladePoints = 0;
    for (int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];

        // Define which way the shaft points to distinguish between
        // upwind and downwind turbines.
        uvShaftDir.append(OverHang[j]/mag(OverHang[j]));

        // Define the vector along the shaft pointing in the
        // direction of the wind.
        uvShaft.append(rotorApex[i] - towerShaftIntersect[i]);
        uvShaft[i] = (uvShaft[i]/mag(uvShaft[i])) * uvShaftDir[i];

        // Define the vector aligned with the tower pointing from
        // the ground to the nacelle.
        uvTower.append(towerShaftIntersect[i] - baseLocation[i]);
        uvTower[i] = uvTower[i]/mag(uvTower[i]);

        // Calculate the width of each actuator section.
        db.append(DynamicList<scalar>(0));
        if (pointDistType[i] == "uniform")
        {
            scalar actuatorWidth = (TipRad[j]-HubRad[j])/numBladePoints[i];
            for (int m = 0; m < numBladePoints[i]; m++)
            {
                db[i].append(actuatorWidth);
            }
        }
        // Add other point distribution types here, such as cosine, tanh.
        else if (pointDistType[i] == "cosine")
        {
            scalar
            anglePerSection =
                Foam::constant::mathematical::pi/numBladePoints[i];
            scalar bladeLength = (TipRad[j]-HubRad[j]);
            for (int m = 0; m < numBladePoints[i]; m++)
            {
                db[i].append
                (
                    -(cos((m + 1)*anglePerSection)
                   - cos(m*anglePerSection))*bladeLength/2.0
                );
            }
        }
        else if (pointDistType[i] == "halfCosine")
        {
            scalar
            anglePerSection =
                0.5*Foam::constant::mathematical::pi/numBladePoints[i];
            scalar bladeLength = (TipRad[j]-HubRad[j]);
            for (int m = 0; m < numBladePoints[i]; m++)
            {
                db[i].append
                (
                    -(
                        cos
                        (
                            0.5*Foam::constant::mathematical::pi
                          + (m+1)*anglePerSection
                        )
                      - cos
                        (
                            0.5*Foam::constant::mathematical::pi
                          + m*anglePerSection
                        )
                     ) * bladeLength
                );
            }
        }

        // Now calculate the actuator section center points for each blade
        // of each turbine in the array.  All blades points will be calculated
        // at zero azimuth (blade pointing up), and then rotated to its correct
        // position before doing a global rotation to the initial azimuth of
        // the rotor.  Also calculate the radius of each point (not including
        // coning).
        bladePoints.append
        (
            List<List<vector> >
            (
                NumBl[j], List<vector>(numBladePoints[i],vector::zero)
            )
        );
        bladeRadius.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
//      Info<<"#BLADE POINTS: "<<endl;  //mark
//      Info<<endl;             //mark
//      Info<<"#Blade 0"<<endl;         //mark
        for (int k = 0; k < NumBl[j]; k++)
        {
            vector root = rotorApex[i];
            scalar beta = PreCone[j][k] - ShftTilt[j];
            root.x() = root.x() + HubRad[j]*Foam::sin(beta);
            root.z() = root.z() + HubRad[j]*Foam::cos(beta);
            scalar dist = 0.0;
            for (int m = 0; m < numBladePoints[i]; m++)
            {
               dist = dist + 0.5*db[i][m];
               bladePoints[i][k][m].x() = root.x() + dist*Foam::sin(beta);
               bladePoints[i][k][m].y() = root.y();
               bladePoints[i][k][m].z() = root.z() + dist*Foam::cos(beta);
               bladeRadius[i][k][m] = HubRad[j] + dist;
               totBladePoints++;
               dist = dist + 0.5*db[i][m];
//             Info<<bladePoints[i][k][m].x()<<" \t";     //mark
//             Info<<bladePoints[i][k][m].y()<<" \t";     //mark
//             Info<<bladePoints[i][k][m].z()<<endl; //mark
            }
            // Apply rotation to get blades, other than blade 1, in the right
            // place.
            if (k > 0)
            {
//              Info<<endl;                   //mark
//              Info<<"#Blade "<<k<<endl;     //mark
                for (int m = 0; m < numBladePoints[i]; m++)
                {
                    bladePoints[i][k][m] =
                        rotatePoint
                        (
                            bladePoints[i][k][m],
                            rotorApex[i],
                            uvShaft[i],
                            (360.0/NumBl[j])*k*degRad
                        );
//              Info<<bladePoints[i][k][m].x()<<" \t";     //mark
//              Info<<bladePoints[i][k][m].y()<<" \t";     //mark
//              Info<<bladePoints[i][k][m].z()<<endl; //mark
                }
            }
        }

        // Define the size of the deltaNacYaw and deltaAzimuth lists and set to
        // zero.
        deltaAzimuth.append(0.0);
        pitch.append(0.0);
        deltaNacYaw.append(0.0);

        // Define the size of the bladeAlignedVectors array and set to zero.
        bladeAlignedVectors.append
        (
            List<List<vector> >(NumBl[j],List<vector>(3,vector::zero))
        );
        chordAlignedVectors.append
        (
            List<List<List<vector> > >
            (
                NumBl[j],
                List<List<vector> >
                (
                    numBladePoints[i],
                    List<vector>(3,vector::zero)
                )
            )
        );
        samplePoints.append
        (
            List<List<List<vector> > >
            (
                NumBl[j],
                List<List<vector> >
                (
                    numBladePoints[i],
                    List<vector>(numSamplePoints[i],vector::zero)
                )
            )
        );

        // Define the size of the samplePointsI array and set to
        // zero.
        samplePointsI.append
        (
            List<List<List<label> > >
            (
                NumBl[j],
                List<List<label> >
                (
                    numBladePoints[i],
                    List<label>(numSamplePoints[i],0)
                )
            )
        );

        // Define the size of the samplePointsU array and set to zero.
        samplePointsU.append
        (
            List<List<List<vector> > >
            (
                NumBl[j],
                List<List<vector> >
                (
                    numBladePoints[i],
                    List<vector>(numSamplePoints[i],vector::zero)
                )
            )
        );

        // Define the size of the circulation sample variables array and set to
        // zero.
        rs.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 0.0))
        );
        qC.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 0.0))
        );
        VmagC.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 0.0))
        );
        alphaC.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 0.0))
        );
        GammaC.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 0.0))
        );
        epsilonC.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 
                epsilon[i]))
        );
        smearingRadiusC.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i], 
                smearRadius[i]))
        );

        // Define the size of the thrust lists and set to zero.
        thrust.append(0.0);
        // Define the size of the aerodynamic torque lists and set to zero.
        torqueRotor.append(0.0);
        // Define the size of the rotor power lists and set to zero.
        powerRotor.append(0.0);

        // Define the size of the coefficient of lift lists and set to zero.
        Cl.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the coefficient of drag lists and set to zero.
        Cd.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the lift/density lists and set to zero.
        lift.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the drag/density lists and set to zero.
        drag.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );

        // Define the size of the bladeForce array and set to zero.
        bladeForce.append
        (
            List<List<vector> >
            (
                NumBl[j],
                List<vector>(numBladePoints[i],vector::zero)
            )
        );
        // Define the size of the bladeForceCorrector array and set to zero.
        bladeForceCorrector.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the axial force/density lists and set to zero.
        axialForce.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the tangential force/density lists and set to zero
        // .
        tangentialForce.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the normal force/density lists and set to zero
        // .
        normalForce.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        // Define the size of the chordwise force/density lists and set to zero
        // .
        chordwiseForce.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );

        // Interpolate to create a global table for co-location points
        // Declare global table for co-location chord and twist
        Chord.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        Twist.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
        TCratio.append
        (
            List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0))
        );
    
        // Interpolate the chord and twist to build the co-location point chord and twist
        for (int k = 0; k < NumBl[j]; k++)
        {
            for (int m = 0; m < numBladePoints[i]; m++)
            {
                // Interpolate the local chord.
                Chord[i][k][m] =
                    interpolate
                    (
                        bladeRadius[i][k][m],
                        DesignStation[j],
                        DesignChord[j]
                    );
                // Interpolate the local twist.
                Twist[i][k][m] =
                    interpolate
                    (
                        bladeRadius[i][k][m],
                        DesignStation[j],
                        DesignTwist[j]
                    );
                // Interpolate the local thickness-to-chord ratio.
                TCratio[i][k][m] =
                    interpolate
                    (
                        bladeRadius[i][k][m],
                        DesignStation[j],
                        DesignTC[j]
                    );

//              if(i==0 and k==0) {Info << Chord[i][k][m] << "\t" << Twist[i][k][m] << endl;}
            }
        }
        
        // Read table of Lift and Drag Coef. and create global table for co-location points
    }

    openOutputFiles();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::actuationLineSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const volVectorField& U = eqn.psi();
    autoPtr<interpolation<vector> >
    Uinterp = interpolation<vector>::New(interpolationDict, U);

    // Update the time step size.
    dt = runTime_.deltaT().value();

    // Update the current simulation time.
    time = runTime_.timeName();

    computeRotSpeed();
    computePitch();
    rotateBlades();
    computeNacYaw();
    yawNacelle();
    computeChordAlignedVectors();
    computeSamplePoints();

    {
        cellFinder *= 0.0;
        clock_t t = ::clock();
        Info<< "sampling velocity ... " << endl;
        forAll(samplePoints, i) // turbine
        {
            forAll(samplePoints[i], j) // blade
            {
                forAll(samplePoints[i][j], k) // actuator point
                {
                    forAll(samplePoints[i][j][k], l) // sample point
                    {
                        label tmpI = -1;
                        tmpI = ms.findCell(samplePoints[i][j][k][l]);
                        vector tmpU = vector::zero;
                        if (tmpI == -1)
                        {
                            tmpI = 0;
                        }
                        else
                        {
                            cellFinder[tmpI] = 1.0;
                            tmpU = U[tmpI];
                            tmpU =
                                Uinterp->interpolate
                                (
                                    samplePoints[i][j][k][l],
                                    tmpI
                                );
                        }
                        label tmpImax = tmpI;
                        reduce(tmpI, sumOp<label>());
                        reduce(tmpImax, maxOp<label>());
                        if (tmpI != tmpImax)
                        {
                            Info<< "WARNING! tmpIsum != tmpImax for samplePoint: " << samplePoints[i][j][k][l] << endl;
                        }
                        else if (tmpI == 0)
                        {
                            Info<< "WARNING! tmpIsum == 0 for samplePoint: " << samplePoints[i][j][k][l] << endl;
                        }
                        reduce(tmpU, sumOp<vector>());
                        samplePointsI[i][j][k][l] = tmpI;
                        if (mag(tmpU) == 0 || tmpI != tmpImax || tmpI == 0)
                        {
                            samplePointsU[i][j][k][l] = rotorUref[i];
                        }
                        else
                        {
                            samplePointsU[i][j][k][l] = tmpU;
                        }
                    }
                }
            }
        }
        t = ::clock() - t;
        Info<< "sampling velocity done (" << (float(t))/CLOCKS_PER_SEC
            << " secs)." << endl;
    }

    computeBladeForceCorrector();
    computeBladeForce();
    computeBodyForce();

    eqn -= -bodyForce;

    // Print turbine output to file.
    outputIndex++;

    if (outputControl == "timeStep")
    {
        if (outputIndex >= outputInterval)
        {
            outputIndex = 0;
            printOutputFiles();
        }
    }
    else if (outputControl == "runTime")
    {
        if ((runTime_.value() - lastOutputTime) >= outputInterval)
        {
            lastOutputTime += outputInterval;
            printOutputFiles();
        }
    }
    else
    {
        printOutputFiles();
    }

    if (debug_)
    {
        printDebug();
    }
}


void Foam::fv::actuationLineSource::writeData(Ostream& os) const
{
}


bool Foam::fv::actuationLineSource::read(const dictionary& dict)
{
    return false;
}


// ************************************************************************* //
