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
#include <string>


void testStdInput(){
    char userinput[100];
    std::cout<<"Input something: ";
    std::cin.getline(userinput,100);
    std::cout<<"The user input '"<<userinput<<"'";
}