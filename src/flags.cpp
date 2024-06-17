/**
 * @file flags.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2022-09-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "flags.h"
bool g_MATRICES_LOADED = false;        //Tells APC whether the picard iteration matrices have already been loaded
bool g_DEBUG_MESSAGES = false;
bool g_KERNELS_LOADED = false;
bool g_DEBUG_PICARD = false;
bool g_DEBUG_INTERPOLATE = false;
bool g_DEBUG_SEGMENTS = false;       //Tells APC whether to store the integrators segment solutions
bool g_BOOTSTRAP_DISABLED = false;   //Tells APC whether to use the bootstrap method
bool g_DEBUG_BOOTSRAP = false;
bool g_DISABLE_REOSCULATION = true; // FIXME: This is a temporary flag to guarantee no reosculation is done. This should be removed once the reosculation is fixed for the bootstrap orbit
bool g_VERBOSE = false;  