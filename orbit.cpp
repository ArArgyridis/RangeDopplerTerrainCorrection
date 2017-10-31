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

#include "orbit.h"

using namespace std;
using namespace boost::posix_time;

PosVector::PosVector():x(0.0), y(0.0), z(0.0) {}

SRGRCoeficientList::SRGRCoeficientList():time(0), groundRangeOrigin(0){}

SRGRCoeficientList::~SRGRCoeficientList(){}

Orbit :: Orbit(orbitVectorType& orbitStateVectors, double &firstLineUTC, double &lineTimeInterval, int& sourceImageHeight) {
    removeRedundantVectors(orbitStateVectors);

    sensorPositionVectors.resize(sourceImageHeight);
    initVelocityTime = orbitStateVectors[0].time;
    dt = (orbitStateVectors.back().time - orbitStateVectors[0].time ) / (orbitStateVectors.size() -1);

    //computing position velocity
    size_t i;
//#pragma omp parallel for private(i) schedule(dynamic, 1500)
    for ( i = 0; i < sourceImageHeight; i++) {
        double tmpTime = firstLineUTC + i*lineTimeInterval;
        getPositionVelocity(tmpTime, sensorPositionVectors[i]);

    }

    numberOfStateVectors = sensorPositionVectors.size();
    //cout << numberOfStateVectors << endl;
}

void Orbit::getPositionVelocity(double &dopplerTime, OrbitData &inputFeature) {
    int nv = 8;
    int i0, iN;
    if ( orbitStateVectors.size() <= nv ) {
        i0 = 0;
        iN = orbitStateVectors.size() - 1;
    }
    else {
        i0 = max( (int) ((dopplerTime - initVelocityTime)/dt - nv/2 + 1), 0);
        iN = min(i0 + nv - 1, (int)(orbitStateVectors.size() - 1));
        i0 = (iN < orbitStateVectors.size() - 1? i0 : iN - nv + 1);
    }

    for (register size_t i = i0; i <=iN; ++i) {
        OrbitData tmpStateVector = orbitStateVectors[i];

        double weight = 1.0;

        for (register size_t j = i0; j <= iN; ++j) {
            if (j != i) {
                double  time2 =  orbitStateVectors[j].time;
                 weight *=  (dopplerTime - time2)*1.0/(tmpStateVector.time - time2);
            }
        }
        // set position
        inputFeature.x_pos += weight*tmpStateVector.x_pos;
        inputFeature.y_pos += weight*tmpStateVector.y_pos;
        inputFeature.z_pos += weight*tmpStateVector.z_pos;

        //set velocity
        inputFeature.x_vel +=weight*tmpStateVector.x_vel;
        inputFeature.y_vel +=weight*tmpStateVector.y_vel;
        inputFeature.z_vel +=weight*tmpStateVector.z_vel;

    }
}

void Orbit ::removeRedundantVectors(orbitVectorType &orbVec) {

    for (register size_t i = 0; i < orbVec.size(); i++)
        if ( i == 0)
            orbitStateVectors.push_back( orbVec[i] );
        else {
            if (orbVec[i].time > orbitStateVectors.back().time  )
                orbitStateVectors.push_back( orbVec[i] );
        }
}


Orbit :: ~Orbit() {}

OrbitData::OrbitData():x_pos(0.0), y_pos(0.0), z_pos(0.0), x_vel(0.0), y_vel(0.0), z_vel(0.0), slant_range(0.0), range_index (0.0), azimuth_index(0.0) {}
