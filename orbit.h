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

#ifndef ORBIT_H
#define ORBIT_H

#include <libxml/parser.h>
#include <libxml/xmlIO.h>
#include <libxml/xinclude.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <ogrsf_frmts.h>
#include <gdal_priv.h>
#include <gdal.h>
#include <ogr_feature.h>
#include <ogr_geometry.h>
#include <algorithm>
#include <iostream>


#include "functions.h"

class SRGRCoeficientList {
public:
    double time, groundRangeOrigin;
    std::vector<double> coeficients;
    SRGRCoeficientList();
    ~SRGRCoeficientList();
};

class PosVector {
public:
    double x, y, z;
    PosVector();
};

class OrbitData {
public:
    double x_pos, y_pos, z_pos, x_vel, y_vel, z_vel, slant_range, range_index, azimuth_index;
    double time;
    OrbitData();
};


typedef std::vector<OrbitData> orbitVectorType;

class RangeDopplerGeocoding;

class Orbit
{
    friend class RangeDopplerGeocoding;
    double dt; //in days
    int numberOfStateVectors;
    orbitVectorType orbitStateVectors, sensorPositionVectors;
    double initVelocityTime;

    void getPositionVelocity(double& dopplerTime, OrbitData& inputFeature);
    void removeRedundantVectors(orbitVectorType& orbVec);

public:
    Orbit(orbitVectorType& orbitStateVectors, double& firstLineUTC, double& lineTimeInterval, int& sourceImageHeight);
    ~Orbit();
};

#endif // ORBIT_H
