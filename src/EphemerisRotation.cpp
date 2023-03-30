/**
 * @file EphemerisRotation.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2023-03-26
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <string>
#include <Orbit.h>
#include "SpiceUsr.h"
#include "EphemerisRotation.h"

void BodyFixedToInertial(vector<double> xF, vector<double> vF, vector<double> &xI, vector<double> &vI, double t, Orbit &orbit)
{
    string fixedCenter = orbit.GetPrimaryBody();
    string InertCenter = orbit.GetInertCenter();
    string FixedFrame = orbit.GetFixedFrame();
    string InertFrame = orbit.GetInertFrame();
    double et = orbit.et(t);
    double FState[6] = {xF[0],xF[1],xF[2],vF[0],vF[1],vF[2]};
    double IState[6];
    double xform[6][6];
    double inertialPrimaryState[6];
    double lt;

    //Get state rotation matrix
    sxform_c(FixedFrame.c_str(),InertFrame.c_str(),et, xform);
    spkezr_c(fixedCenter.c_str(),et,InertFrame.c_str(),"LT+S",InertCenter.c_str(),inertialPrimaryState,&lt);
    mxvg_c(xform,FState,6,6,IState);
    for(int i=0;i<6;i++){
        IState[i]+=inertialPrimaryState[i];
    }
    xI.insert(xI.begin(),IState,IState+3);
    vI.insert(vI.begin(),IState+3,IState+6);
    
}
void BodyFixedToInertial(double *xF, double *vF, double *xI, double *vI, double t, Orbit &orbit)
{
    string fixedCenter = orbit.GetPrimaryBody();
    string InertCenter = orbit.GetInertCenter();
    string FixedFrame = orbit.GetFixedFrame();
    string InertFrame = orbit.GetInertFrame();
    double et = orbit.et(t);
    double FState[6] = {xF[0],xF[1],xF[2],vF[0],vF[1],vF[2]};
    double IState[6];
    double xform[6][6];
    double inertialPrimaryState[6];
    double lt;

    //Get state rotation matrix
    sxform_c(FixedFrame.c_str(),InertFrame.c_str(),et, xform);
    spkezr_c(fixedCenter.c_str(),et,InertFrame.c_str(),"LT+S",InertCenter.c_str(),inertialPrimaryState,&lt);
    mxvg_c(xform,FState,6,6,IState);
    //Add 
    for(int i=0;i<3;i++){
        xI[i]=IState[i]+inertialPrimaryState[i];
        vI[i]=IState[i+3]+inertialPrimaryState[i+3];
    }
}
void BodyFixedAccelerationToInertial(double *aF, double *aI, double t, Orbit &orbit)
{
    string fixedCenter = orbit.GetPrimaryBody();
    string InertCenter = orbit.GetInertCenter();
    string FixedFrame = orbit.GetFixedFrame();
    string InertFrame = orbit.GetInertFrame();
    double et = orbit.et(t);
    double FState[3] = {aF[0],aF[1],aF[2]};
    double IState[3];
    double xform[3][3];
    double lt;

    //Get state rotation matrix
    pxform_c(FixedFrame.c_str(),InertFrame.c_str(),et, xform);
    mxvg_c(xform,FState,3,3,IState);
    for(int i=0;i<3;i++){
        aI[i]=IState[i];
    }
}
void InertialToBodyFixed(vector<double> xI, vector<double> vI, vector<double> &xF, vector<double> &vF, double t, Orbit &orbit)
{
    string fixedCenter = orbit.GetPrimaryBody();
    string InertCenter = orbit.GetInertCenter();
    string FixedFrame = orbit.GetFixedFrame();
    string InertFrame = orbit.GetInertFrame();
    double et = orbit.et(t);
    double IState[6] = {xI[0],xI[1],xI[2],vI[0],vI[1],vI[2]};
    double FState[6];
    double xform[6][6];
    double inertialPrimaryState[6];
    double lt;

    //Substract moon state from inertial satellite state to get state relative to the moon
    spkezr_c(fixedCenter.c_str(),et,InertFrame.c_str(),"LT+S",InertCenter.c_str(),inertialPrimaryState,&lt);
    for(int i=0;i<6;i++){
    IState[i]-=inertialPrimaryState[i];
    }
    //Get state rotation matrix and rotate from inertial frame to fixed frame.
    sxform_c(InertFrame.c_str(),FixedFrame.c_str(),et, xform);
    mxvg_c(xform,IState,6,6,FState);
    //copy state to position and velocity vector
    xF.insert(xF.begin(),FState,FState+3);
    vF.insert(vF.begin(),FState+3,FState+6);
    
}

void InertialToBodyFixed(double *xI, double *vI, double *xF, double *vF, double t, Orbit &orbit)
{
    string fixedCenter = orbit.GetPrimaryBody();
    string InertCenter = orbit.GetInertCenter();
    string FixedFrame = orbit.GetFixedFrame();
    string InertFrame = orbit.GetInertFrame();
    double et = orbit.et(t);
    double IState[6] = {xI[0],xI[1],xI[2],vI[0],vI[1],vI[2]};
    double FState[6];
    double xform[6][6];
    double inertialPrimaryState[6];
    double lt;

    //Substract moon state from inertial satellite state to get state relative to the moon
    spkezr_c(fixedCenter.c_str(),et,InertFrame.c_str(),"LT+S",InertCenter.c_str(),inertialPrimaryState,&lt);
    for(int i=0;i<6;i++){
    IState[i]-=inertialPrimaryState[i];
    }
    //Get state rotation matrix and rotate from inertial frame to fixed frame.
    sxform_c(InertFrame.c_str(),FixedFrame.c_str(),et, xform);
    mxvg_c(xform,IState,6,6,FState);
    //copy state to position and velocity vector
    for(int i=0;i<3;i++){
        xF[i]=FState[i];
        vF[i]=FState[i+3];
    }
}
