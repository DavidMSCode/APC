/**
 * @file preferences.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief Functions for storing user preferences for APC
 * @version 0.1
 * @date 2022-08-24
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <iostream>
#include <fstream>
#include <cstdint>
#include <filesystem>
#include <string>
#include <cstdlib>

#include "preferences.h"
namespace fs = std::filesystem;

void testStdInput(){
    char userinput[100];
    std::cout<<"Input something: ";
    std::cin.getline(userinput,100);
    std::cout<<"The user input '"<<userinput<<"'";
}

bool userPreferenceExists(){
    bool exists;
    //Get preference file location assuming 
    const fs::path  preference_path = preferenceFilePath();
    if(fs::exists(preference_path)){
        exists = true;
    }
    else{
        exists = false;
    }
    return exists;
}

std::string preferenceFilePath(){
    //get path to $HOME defined by environment variable
    std::string home_path = std::string(std::getenv("HOME"));
    const char *preference_path;
    //Check for operating system
    #if (_WIN32 || _WIN64 || __CYGWIN__)
        //Windows
        preference_path = &WINDOWS_PREF_LOCATION[0];
    #elif __MACH__
        //Mac OS
        preference_path = &MACOS_PREF_LOCATION[0];
    #elif __linux__
        //some kind of Linux
        preference_path = &LINUX_PREF_LOCATION[0];
    #else
        //None of the above assume default location;
        preference_path = &DEFAULT_PREF_LOCATION[0];
    #endif
    //combine home path and preference file path to get full path to file
    std::string full_path = home_path +"/"+ std::string(preference_path);
    // std::cout<<full_path+"endl";
    return full_path;
}
