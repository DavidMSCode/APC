/**
 * @file preferences.h
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2022-08-24
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __PREF__
#define __PREF__

//Default Values
const char* const DEFAULT_KERNELS = "de440.bsp"; //Spice kernel to be loaded
const char* const DEFAULT_TLS = "naif0012.tls"; //default leap second kernel
//Relative paths (from $HOME) to preference files for different OSes
const char* const MACOS_PREF_LOCATION = "Library/Preferences/APC.xml"; //mac os user preference file location
const char* const LINUX_PREF_LOCATION = ".APC.xml";    //Linux user preference file location
const char* const WINDOWS_PREF_LOCATION = ".APC.xml"; //WIndows user preference file location
const char* const DEFAULT_PREF_LOCATION = ".APC.xml"; //default location when os is not known

/**
 * @brief test pybind11 interaction with std::cin
 * 
 */
void testStdInput();

/**
 * @brief returns true if user preference file exists
 * 
 */
bool userPreferenceExists();

/**
 * @brief Returns the file path location for the current OS (Mac, Windows, Linux)
 * 
 * @return const char* a pointer to a const character array defining the preferenceFileLocation
 */
std::string preferenceFilePath();


#endif