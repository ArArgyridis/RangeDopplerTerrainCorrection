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

#include "functions.h"

using namespace boost::posix_time;
using namespace std;

double computePolynomialValue(double &x, std::vector<double>& coefs) {
    double v = 0.0;
    for (int i = (int)coefs.size() - 1; i > 0; i--)
        v = (v + coefs[i])*x;
    return v + coefs[0];
}

bool createDirectory(std::string &directory) {
    struct stat info;
    if( stat( directory.c_str(), &info ) != 0 ) {
        cout << "creating output directory: " << directory << endl;
        size_t pre=0, pos;
        int mdret;

        if(directory.back() != '/')
            directory+='/';

        bool checkDir = false;
        while( (pos=directory.find_first_of('/',pre)) != std::string::npos) {
                string tmpDir = directory.substr(0, pos++);
                pre = pos;
                checkDir = mkdir(tmpDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
        }

        if (checkDir == 0) {
            cout << "output directory created!\n";
            return true;
        }
        else{
            cerr <<"error: unable to create output directory\n";
            return false;
        }
    }
    else if( info.st_mode & S_IFDIR ) {
        cout << "Warning: output directory already exists!\n";
        return true;
    }
}

void splitString(const string &inStr, char &splitChar, vector<string> &outVec) {
    string seg;
    stringstream tmpInStr(inStr);
    while(std::getline(tmpInStr, seg, splitChar))
    {
        outVec.push_back(seg);
    }
}

double stringToDate(string &inDateTime) {
    time_input_facet* facet = new time_input_facet("%d-%b-%Y %H:%M:%s");
    stringstream ss(inDateTime);
    ss.imbue(locale(locale(), facet));
    ptime pt;
    ss >> pt;
    //cout <<setprecision(20) << (long long)((pt - ptime { {1970,1,1}, {} }).total_nanoseconds()) << endl;
    return (pt - ptime { {1970,1,1}, {} }).total_nanoseconds();

}
