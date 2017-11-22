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

#ifndef RANGEDOPPLERGEOCODING_H
#define RANGEDOPPLERGEOCODING_H

#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <boost/date_time.hpp>
#include <sys/types.h>
#include "orbit.h"

const double nullDopplerTime = -9999.0;
const double lightSpeedInMetersPerDay = 299792458.0 * 86400.0;

class RangeDopplerGeocoding
{
    //private variables
    bool skipBistaticCorrection, srgrFlag, isPolsar, nearRangeOnLeft  ;
    int numOfInputBands, sourceImageWidth, sourceImageHeight, targetImageWidth, targetImageHeight, demCols, demRows, sourceSRSEPSG, targetSRSEPSG, numOfSRGRCoeficientVectors;
    double waveLength, nearEdgeSlantRange, avgSceneHeight, xMin, xMax, yMin, yMax, rangeSpacing, azimuthSpacing, pixelSpacing, pixelSpacingInMeters, pixelSpacingInDegrees, semiMajorAxis, *demGeoTransform;
    orbitVectorType orbitVector; //vector containing the state vectors
    std::string  fileName, inDEMFileName, inPath, mission, outputFile;
    std::map<std::string, float*> inputBands;
    std::map<std::string, std::string> outputImageMap;
    std::vector<SRGRCoeficientList> SRGRCoeficientListVector;
    double firstLineTime, lastLineTime, lineTimeInterval;
    Orbit *orbit;
    OGRSpatialReference sourceSRS, targetSRS;
    int xBlockSize, yBlockSize;

    //private functions
    double computeGroundRange(double& groundRangeOrigin, double& slantRange, std::vector<double>& coeficients);
    double computeSlantRange(double& zeroDopplerTime, double& x, double&y, double& z, OrbitData* inSensorPositionFeature );
    double computeRangeIndex(double& zeroDopplerTime, double& slantRange);
    void computeTargetExtents();
    void createOutputDataset(int& xBlockCount, int& yBlockCount); //creates new dataset and returns the number of blocks
    double getDopplerFrequency(double& x, double&y, double& z, int& id );
    double getDoubleProperty(xmlChar* xpath);
    double getHeight(double& x, double& y, int& bufferWidth, int& colOffset, int& rowOffset, float* inDEMDataBuffer);
    double getEarthPointZeroDopplerTime(double& x, double&y, double& z);
    double getPixelValue(OrbitData* sensorPositionFeature, const std::string& bandName);
    std::string getStringProperty(xmlChar* xpath);
    int getIntProperty(xmlChar* xpath);
    bool getPosition(double&x, double& y, double& z, OrbitData* sensorPositionFeature);

    bool isValidCell(OGRFeature* sensorPositionFeature);
    void readBands();


    xmlXPathObjectPtr getNodeSet (xmlChar* xpath); //function to get XPath results from metadata doc
    xmlDocPtr metaDataFile; //metadata file having the information

    void parseMetaDataFile();
public:
    void getTiePointGrid();
    void resampleImage(float& pixelSize);
    RangeDopplerGeocoding(char* inFile, char* inDEM, int targetsrsCode, char* outFile);
    ~RangeDopplerGeocoding();
    void setTargetExtents(double &xmin, double& ymin, double& xmax, double &ymax);
};

#endif // RANGEDOPPLERGEOCODING_H
