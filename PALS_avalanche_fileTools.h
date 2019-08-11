#ifndef PALS_avalanche_fileTools_H
#define PALS_avalanche_fileTools_H

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <cmath>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <iomanip>

class FileTools
{
	public:
		FileTools();
		FileTools( std::string fitDetailsFile );
		
		std::string NumberToChar( double Number, int Precision);    //Returning number as a string with Precision given by Precision argument -> Precision = 2, then Number returned up to 2 decimals
		bool PathCheck( std::string Path );                         //Checking if the Directory given by Path argument exists
		bool FileCheck(const std::string& NameOfFile);              //Checking if the File given by NameOfFile exists
		bool CheckTheLine( char FirstCharacter, char TestCharacter, unsigned NumberOfLine );    //Checking if the line number NumberOfLine starts with TestCharacter
		std::vector<std::string> GetFitDetails();
		std::vector<int> GetComponentsMultiplicities( std::string SizesToExtract );
		std::vector<std::string> GetComponentDetails( std::string DetailsToExtract );
        
	private:
		std::string FitDetailsFile;
};

#endif
