/**
 * @file flags.h
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2022-08-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __FLAGS__
#define __FLAGS__
    extern bool g_MATRICES_LOADED;        //Tells APC whether the picard iteration matrices have already been loaded
    extern bool g_DEBUG_MESSAGES;
    extern bool g_KERNELS_LOADED;       //Tells APC whether the spice kernels are currently loaded
    extern bool g_DEBUG_PICARD;
    extern bool g_DEBUG_INTERPOLATE;
    extern bool g_DEBUG_SEGMENTS;       //Tells aAPC whether to store the integrators segment solutions
    extern bool g_BOOTSTRAP_DISABLED;   //Tells APC whether to use the bootstrap method
    extern bool g_DEBUG_BOOTSRAP;
    extern bool g_DISABLE_REOSCULATION;
#endif