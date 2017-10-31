/*This is part of the program RangeDoppler Geocoding
   Copyright (C) 2017  Argyros Argyridis arargyridis@gmail.com

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <boost/algorithm/string.hpp>
#include <boost/date_time/gregorian/gregorian_io.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <istream>
#include <sys/stat.h>
#include <sys/dir.h>

bool createDirectory(std::string &directory);
double computePolynomialValue(double &x, std::vector<double>& coefs);
void splitString(const std::string& inStr, char& splitChar, std::vector<std::string>& outVec);
double stringToDate(std::string& inDateTime);
#endif // FUNCTIONS_H
