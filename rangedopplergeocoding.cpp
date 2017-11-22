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

#include "rangedopplergeocoding.h"


using namespace std;
using namespace boost::gregorian;
using namespace  boost::posix_time;



double RangeDopplerGeocoding::computeGroundRange(double &groundRangeOrigin, double &slantRange, std::vector<double> &coeficients) {

    double lowerBound = groundRangeOrigin;
    double lowerBoundSlantRange = computePolynomialValue(lowerBound, coeficients);

    if (slantRange < lowerBoundSlantRange)
        return -1.0;

    double upperBound = groundRangeOrigin + sourceImageWidth*rangeSpacing;
    double upperBoundSlantRange =  computePolynomialValue(upperBound, coeficients);

    if (slantRange > upperBoundSlantRange)
        return -1.0;

    double midSlantRange;

    while (upperBound - lowerBound > 0) {
        double mid = (lowerBound + upperBound) / 2.0;
        midSlantRange = computePolynomialValue(mid, coeficients);
        double a = midSlantRange - slantRange;

        if ((a > 0 && a < 0.1) || (a <= 0.0 && 0.0 - a < 0.1)) {
            return mid;
        } else if (midSlantRange < slantRange) {
            lowerBound = mid;
        } else if (midSlantRange > slantRange) {
            upperBound = mid;
        }
    }
    return -1.0;
}

double RangeDopplerGeocoding::computeRangeIndex(double &zeroDopplerTime, double &slantRange) {

    if (zeroDopplerTime < min(firstLineTime, lastLineTime) || zeroDopplerTime > max(firstLineTime, lastLineTime))
        return -1.0;

    if (srgrFlag) {
        double groundRange;
        if( SRGRCoeficientListVector.size() == 1 ) {
            int idx = 0;
            groundRange = computeGroundRange(SRGRCoeficientListVector[0].groundRangeOrigin, slantRange, SRGRCoeficientListVector[0].coeficients);
            if (groundRange < 0.0)
                return -1.0;
            else
                return (groundRange - SRGRCoeficientListVector[0].groundRangeOrigin) / rangeSpacing;
        }
        else {

            int idx = 0;
            for (int i = 0; i < SRGRCoeficientListVector.size() && zeroDopplerTime >= SRGRCoeficientListVector[i].time; i++)
                idx = i;


            vector<double> tmpCoefs;
            tmpCoefs.resize( SRGRCoeficientListVector[idx].coeficients.size() );
            if (idx == numOfSRGRCoeficientVectors -1)
                idx--;

            double mu = (zeroDopplerTime - SRGRCoeficientListVector[idx].time) / ( SRGRCoeficientListVector[idx+1].time - SRGRCoeficientListVector[idx].time );

            for (register size_t i = 0; i < SRGRCoeficientListVector[idx].coeficients.size(); i++) {
                tmpCoefs[i] = (1- mu ) * SRGRCoeficientListVector[idx].coeficients[i] + mu*SRGRCoeficientListVector[idx+1].coeficients[i];
            }

            groundRange = computeGroundRange(SRGRCoeficientListVector[idx].groundRangeOrigin, slantRange, tmpCoefs);

            if (groundRange < 0.0)
                return -1.0;
            else
                return(groundRange - SRGRCoeficientListVector[idx].groundRangeOrigin) / rangeSpacing;
        }
    }
    else { //slant range image
        return (slantRange - nearEdgeSlantRange) / rangeSpacing;
    }
}

double RangeDopplerGeocoding::computeSlantRange(double &zeroDopplerTime, double &x, double &y, double &z, OrbitData *inSensorPositionFeature ) {

    orbit->getPositionVelocity(zeroDopplerTime, *inSensorPositionFeature);
    double xDiff = inSensorPositionFeature->x_pos - x;
    double yDiff = inSensorPositionFeature->y_pos - y;
    double zDiff = inSensorPositionFeature->z_pos -  z;

    return sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);
}

void RangeDopplerGeocoding::computeTargetExtents() {
    OGRCoordinateTransformation *transfTargetSource;
    transfTargetSource = OGRCreateCoordinateTransformation(&sourceSRS, &targetSRS);

    transfTargetSource->Transform(1, &xMin, &yMax);
    transfTargetSource->Transform(1, &xMax, &yMin);
}

void RangeDopplerGeocoding::createOutputDataset(int &xBlockCount, int &yBlockCount) {
    //creating output directory

    string driver = "GTiff";

    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");

    double *outGeoTransform;
    outGeoTransform = new double[6];
    outGeoTransform[0] = xMin;
    outGeoTransform[1] = pixelSpacing;
    outGeoTransform[2] = 0.0;
    outGeoTransform[3] = yMax;
    outGeoTransform[4] = 0.0;
    outGeoTransform[5] = -pixelSpacing;

    char *outProjection = nullptr;
    targetSRS.exportToWkt(&outProjection);

    //checking if path is correct
    char splitChar = '/';
    vector<string> splitOutFile;
    splitString(outputFile, splitChar, splitOutFile);

    string outputPath = "";
    for(size_t i = 0; i < splitOutFile.size()  - 1; i++)
        outputPath += splitOutFile[i] + "/";

    bool checkDir = false;
    checkDir = createDirectory(outputPath);
    if (!checkDir)
            exit (1);

    cout <<"Creating output dataset: " << outputFile << endl;
    GDALDataset *tmpDataset;
    tmpDataset = poDriver->Create(outputFile.c_str(), targetImageWidth, targetImageHeight, inputBands.size(), GDT_Float32, NULL);
    tmpDataset->SetGeoTransform(outGeoTransform);
    tmpDataset->SetProjection(outProjection);

    //determining block size
    if (xBlockSize == 0) {
        tmpDataset->GetRasterBand(1)->GetBlockSize(&xBlockSize, &yBlockSize);
        xBlockCount  = (targetImageWidth + (xBlockSize - 1))/xBlockSize;
        yBlockCount  = (targetImageHeight + (yBlockSize -1 ))/yBlockSize;
    }

    int i = 1;
    for (auto& band : inputBands) {
        //check if output directory exists
        vector<string> splitName;
        char dot = '.';
        splitString(band.first, dot, splitName);
        tmpDataset->GetRasterBand(i++)->SetDescription( splitName[0].c_str() );
    }

    GDALClose(tmpDataset);

    delete[] outGeoTransform;
    outGeoTransform = nullptr;
    CPLFree(outProjection);

}

double RangeDopplerGeocoding::getDopplerFrequency(double &x, double &y, double &z, int &id) {

    OrbitData tmpFt = orbit->sensorPositionVectors[id];

    double xDiff = x - tmpFt.x_pos;
    double yDiff = y - tmpFt.y_pos;
    double zDiff = z - tmpFt.z_pos;
    double distance = sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);

    return 2.0 * ( tmpFt.x_vel*xDiff + tmpFt.y_vel*yDiff + tmpFt.z_vel*zDiff  ) / (distance * waveLength);

}

double RangeDopplerGeocoding::getHeight( double &x, double &y, int &bufferWidth, int &colOffset, int &rowOffset, float* inDEMDataBuffer) {
    //converting geodesic to pixel coordinates
    int col = int(-(demGeoTransform[2]*(demGeoTransform[3] - y) - demGeoTransform[0]*demGeoTransform[5] + demGeoTransform[5]*x)/(demGeoTransform[2]*demGeoTransform[4] - demGeoTransform[1]*demGeoTransform[5])) - colOffset;
    int row = int( (demGeoTransform[1]*(demGeoTransform[3] - y) - demGeoTransform[0]*demGeoTransform[4] + demGeoTransform[4]*x)/(demGeoTransform[2]*demGeoTransform[4] - demGeoTransform[1]*demGeoTransform[5])) - rowOffset;

    return inDEMDataBuffer[ bufferWidth*row +col ];
}

    xmlXPathObjectPtr RangeDopplerGeocoding::getNodeSet (xmlChar *xpath) {

    xmlXPathContextPtr context;
    xmlXPathObjectPtr result;

    context = xmlXPathNewContext(metaDataFile);
    if (context == NULL) {
        printf("Error in xmlXPathNewContext\n");
        return NULL;
    }

    result = xmlXPathEvalExpression(xpath, context);
    xmlXPathFreeContext(context);
    if (result == NULL) {
        printf("Error in xmlXPathEvalExpression\n");
        return NULL;
    }

    if(xmlXPathNodeSetIsEmpty(result->nodesetval)){
        xmlXPathFreeObject(result);
        printf("No result\n");
        return NULL;
    }

    return result;
}

double RangeDopplerGeocoding::getDoubleProperty(xmlChar* xpath) {
    xmlXPathObjectPtr result = getNodeSet(xpath);
    return  stold(reinterpret_cast <char*>(xmlNodeGetContent(result->nodesetval->nodeTab[0]->xmlChildrenNode) ) );
}

int RangeDopplerGeocoding::getIntProperty(xmlChar* xpath) {
    xmlXPathObjectPtr result = getNodeSet(xpath);
    return  stoi(reinterpret_cast <char*>(xmlNodeGetContent(result->nodesetval->nodeTab[0]->xmlChildrenNode) ) );
}

string RangeDopplerGeocoding::getStringProperty(xmlChar* xpath) {
    xmlXPathObjectPtr result = getNodeSet(xpath);
    return  string( reinterpret_cast <char*>( xmlNodeGetContent(result->nodesetval->nodeTab[0]->xmlChildrenNode) ) );
}

double RangeDopplerGeocoding::getPixelValue(OrbitData *sensorPositionFeature, const string &bandName) {
    return inputBands[bandName][ sourceImageWidth*(int)sensorPositionFeature->azimuth_index+(int)sensorPositionFeature->range_index];
}

void RangeDopplerGeocoding::getTiePointGrid() {

    double yFirstNear    =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='first_near_lat']");
    double xFirstNear    =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='first_near_long']");

    double yFirstFar       =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='first_far_lat']");
    double xFirstFar       =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='first_far_long']");

    double yLastNear    =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='last_near_lat']");
    double xLastNear    =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='last_near_long']");

    double yLastFar       =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='last_far_lat']");
    double xLastFar       =  getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='last_far_long']");

    xMin = min(min(min(xFirstNear, xLastNear), xFirstFar), xLastFar);
    xMax = max(max(max(xFirstNear, xLastNear), xFirstFar), xLastFar);

    yMin = min(min(min(yFirstNear, yLastNear), yFirstFar), yLastFar);
    yMax = max(max(max(yFirstNear, yLastNear), yFirstFar), yLastFar);

    //computing new output grid
    computeTargetExtents();
}

void RangeDopplerGeocoding::parseMetaDataFile() {

    //retrieving mission type
    mission = getStringProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='MISSION']");
    boost::to_upper(mission);
    if (mission.find("SENTINEL") != string::npos || mission.find("TSX") != string::npos || (mission == "RS2") || (mission.find("CSKS") != string::npos ) )
        skipBistaticCorrection = true;

    //retrieving srgr_flag
    srgrFlag = getIntProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='srgr_flag']");

    if (srgrFlag) { //GET SRGR_COEFICIENTS
        xmlXPathObjectPtr result = getNodeSet((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDElem[@name='SRGR_Coefficients']");
        xmlNodeSetPtr nodeset = result->nodesetval;
        for (xmlNodePtr node =  xmlNextElementSibling(nodeset->nodeTab[0]->xmlChildrenNode); node != nullptr; node =  xmlNextElementSibling(node)) {
            SRGRCoeficientList tmpList;
            for (xmlNodePtr subElement = xmlFirstElementChild(node); subElement != nullptr; subElement  = xmlNextElementSibling(subElement) ) {
                if ( xmlStrEqual(xmlGetProp(subElement, (xmlChar*)"name"), (xmlChar*)"zero_doppler_time") ) {
                    string tmpDate = string(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
                    tmpList.time = stringToDate(tmpDate);
                }
                else if ( xmlStrEqual(xmlGetProp(subElement, (xmlChar*)"name"), (xmlChar*)"ground_range_origin"))
                    tmpList.groundRangeOrigin = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
                else if(  string((char*)xmlGetProp(subElement, (xmlChar*)"name")).find("coefficient" ) != string::npos  ) {
                    tmpList.coeficients.push_back( stod(reinterpret_cast<char*>(xmlNodeGetContent(xmlFirstElementChild(subElement)))));
                }
            }
            SRGRCoeficientListVector.push_back(tmpList);
        }
        numOfSRGRCoeficientVectors = SRGRCoeficientListVector.size();
    }
    else
        nearEdgeSlantRange = getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='slant_range_to_first_pixel']");

    //retrieving radar frequency
    waveLength =getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='radar_frequency']" );

    //retrieving range spacing
    rangeSpacing = getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='range_spacing']");

    //retrieving first_line_time and last_line_time
    string tmpDate = getStringProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='first_line_time']");
    firstLineTime = stringToDate(tmpDate);
    tmpDate = getStringProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='last_line_time']");
    lastLineTime = stringToDate(tmpDate);

    //retrieving average scene height - used for retro-calibration or when useAvgSceneHeight is true
    avgSceneHeight = getDoubleProperty((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='avg_scene_height']");

    //retrieving state vectors from metadata file
    xmlXPathObjectPtr result = getNodeSet((xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDElem[@name='Orbit_State_Vectors']");
    xmlNodeSetPtr nodeset = result->nodesetval;
    for (xmlNodePtr node =  xmlNextElementSibling(nodeset->nodeTab[0]->xmlChildrenNode); node != nullptr; node =  xmlNextElementSibling(node)) {
        OrbitData tmpOrbit;
        for (xmlNodePtr subElement = xmlFirstElementChild( node ); subElement != nullptr; subElement =  xmlNextElementSibling(subElement) ) {
            string tmpProp = string(reinterpret_cast<char*>(xmlGetProp(subElement, (xmlChar*) "name")));
            if ( tmpProp == "x_pos" )
                tmpOrbit.x_pos = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
            else if ( tmpProp == "y_pos" )
                tmpOrbit.y_pos = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
            else if ( tmpProp == "z_pos" )
                tmpOrbit.z_pos = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
            else if ( tmpProp == "x_vel" )
                tmpOrbit.x_vel = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
            else if ( tmpProp == "y_vel" )
                tmpOrbit.y_vel = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
            else if ( tmpProp == "z_vel" )
                tmpOrbit.z_vel = stod(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
            else if ( tmpProp == "time" ) {
                string tmpDate = string(reinterpret_cast<char*>(xmlNodeGetContent(subElement)));
                tmpOrbit.time = stringToDate(tmpDate);
            }
        }
        orbitVector.push_back(tmpOrbit);
    }

    //retrieving image width and height
    sourceImageWidth = getIntProperty( (xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='num_samples_per_line']" );
    sourceImageHeight = getIntProperty( (xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='num_output_lines']" );
    isPolsar = getIntProperty( (xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='polsar_data']" );

    lineTimeInterval = (lastLineTime - firstLineTime) / (sourceImageHeight - 1);
    //getting azimuth and range spacing
    azimuthSpacing = getDoubleProperty( (xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='azimuth_spacing']" );
    rangeSpacing = getDoubleProperty( (xmlChar*)"//Dataset_Sources/MDElem[@name='metadata']/MDElem[@name='Abstracted_Metadata']/MDATTR[@name='range_spacing']" );

    if (!srgrFlag) {
        GDALDataset *incidentAngle;
        string incidentAngleFile = inPath + "tie_point_grids/incident_angle.img";
        incidentAngle = (GDALDataset*)GDALOpen( incidentAngleFile.c_str() , GA_ReadOnly );
        float *incidentAngleCenter;
        incidentAngleCenter = new float;
        int xc, yc;
        xc = incidentAngle->GetRasterXSize()/2;
        yc = incidentAngle->GetRasterYSize()/2;
        incidentAngle->GetRasterBand(1)->RasterIO( GF_Read, xc, yc, 1, 1, incidentAngleCenter, 1, 1, GDT_Float32, 0, 0);
        pixelSpacingInMeters = rangeSpacing/sin(*incidentAngleCenter*M_PI/180.0);
        delete incidentAngleCenter;
        GDALClose(incidentAngle);
    }
    else { //check here what is happening for pixel spacing

    }
    pixelSpacingInMeters = max(azimuthSpacing, pixelSpacingInMeters);
    //getting elispoid semi-major axis
    pixelSpacingInDegrees = pixelSpacingInMeters / 6378137 * 180.0 /M_PI;

    //retrieving band names
    result = getNodeSet((xmlChar*)"//Data_Access");
    nodeset = result->nodesetval;
    //getting input path
    vector<string> splitPath;
    char delimiter = '.';
    splitString(fileName, delimiter, splitPath);

    for (xmlNodePtr node =  xmlNextElementSibling(nodeset->nodeTab[0]->xmlChildrenNode); node != nullptr; node =  xmlNextElementSibling(node)) {
        if ( xmlStrEqual( node->name, (xmlChar*)"Data_File") ) {
            string inBand = (reinterpret_cast<char*>(xmlGetProp(xmlFirstElementChild(node), (xmlChar*)"href")));
            inBand.replace(inBand.length()-3, inBand.length()-1, "img");
            vector<string> splitBand, splitName;
            delimiter = '/';
            splitString(inBand, delimiter, splitBand);
            delimiter = '.';
            splitString(splitBand.back(), delimiter, splitName);
            inputBands[splitName.front() ] = nullptr;
        }
    }
}

void RangeDopplerGeocoding::readBands() {
    for(auto& band : inputBands) {
        cout << "Reading" <<" " << inPath + band.first + ".img" << endl;
        GDALDataset *tmpDataset;
        tmpDataset = (GDALDataset*)( GDALOpen((inPath + band.first + ".img").c_str(), GA_ReadOnly) );
        band.second = new float[sourceImageHeight*sourceImageWidth];
        tmpDataset->GetRasterBand(1)->AdviseRead(0, 0, sourceImageWidth, sourceImageHeight, sourceImageWidth, sourceImageHeight, GDT_Float32, NULL);
        tmpDataset->GetRasterBand(1)->RasterIO( GF_Read, 0, 0, sourceImageWidth, sourceImageHeight, band.second, sourceImageWidth, sourceImageHeight, GDT_Float32, 0, 0);
        GDALClose(tmpDataset);
    }
}

void RangeDopplerGeocoding::resampleImage(float &pixelSize) {
    //checking if the pixel size is to be determined automatically or not
    if (pixelSize != -1.0)
        pixelSpacing = pixelSize;
    else
        targetSRS.IsProjected() ? pixelSpacing = pixelSpacingInMeters : pixelSpacing = pixelSpacingInDegrees;

    targetImageWidth = round((xMax - xMin)/pixelSpacing);
    targetImageHeight = round((yMax - yMin)/pixelSpacing);

    //reading DEM
    GDALDataset *inDEMData;
    inDEMData = (GDALDataset*)( GDALOpen(inDEMFileName.c_str(), GA_ReadOnly) );

    demCols = inDEMData->GetRasterXSize();
    demRows = inDEMData->GetRasterYSize();
    demGeoTransform = new double[6];
    inDEMData->GetGeoTransform(demGeoTransform);

    GDALClose(inDEMData);

    int xBlockCount, yBlockCount;
    createOutputDataset(xBlockCount, yBlockCount);

    cout << "Starting resampling....\n";

    #pragma omp parallel
    {
        OGRSpatialReference geocentricSRS;
        geocentricSRS.importFromEPSG(4978);
        OGRCoordinateTransformation *transfTargetGeocentric, *tmpTransf;
        transfTargetGeocentric = OGRCreateCoordinateTransformation(&targetSRS, &geocentricSRS);
        tmpTransf = OGRCreateCoordinateTransformation(&targetSRS, &sourceSRS);

        //creating a dataset for each thread
        GDALDataset *inDEMData2, *outputDataset;
        inDEMData2 = (GDALDataset*)( GDALOpen(inDEMFileName.c_str(), GA_ReadOnly) );
        outputDataset = (GDALDataset*)( GDALOpen(outputFile.c_str(), GA_Update) );

        map<string, GDALRasterBand*> outBands;
        for (register int i = 1; i < inputBands.size() + 1; i++) {
            GDALRasterBand *tmpBand;
            tmpBand = outputDataset->GetRasterBand(i);
            outBands[ string(tmpBand->GetDescription()) ] = tmpBand;
        }


        GDALRasterBand *inDEMBand;
        inDEMBand = inDEMData2->GetRasterBand(1);

        double x, y;
        map<string, float*> bandBuffer;
        int xBlock,  yBlock;

        for ( auto& band : inputBands)
            bandBuffer[band.first] = new float[xBlockSize*yBlockSize];


        int previousXBlockSize = xBlockSize, previousYBlockSize = yBlockSize;
        #pragma omp for schedule(dynamic, 1) collapse(2)
        for (xBlock = 0; xBlock < xBlockCount; xBlock++) { //for each xblock
            for (yBlock = 0; yBlock < yBlockCount; yBlock++) { // for each yblock

                //getting actual block size
                int actualXBlockSize, actualYBlockSize;
                outputDataset->GetRasterBand(1)->GetActualBlockSize(xBlock, yBlock, &actualXBlockSize, &actualYBlockSize);

                for ( auto& band : inputBands) {
                    if ( (previousXBlockSize != actualXBlockSize) || (previousYBlockSize != actualYBlockSize) ) {
                        delete[] bandBuffer[band.first];
                        bandBuffer[band.first] = nullptr;
                        bandBuffer[band.first] = new float[actualXBlockSize*actualYBlockSize];
                        fill(bandBuffer[band.first], bandBuffer[band.first] + actualXBlockSize*actualYBlockSize, 0.0);
                        previousXBlockSize = actualXBlockSize;
                        previousYBlockSize = actualYBlockSize;
                    }
                    else {
                        fill(bandBuffer[band.first], bandBuffer[band.first] + actualXBlockSize*actualYBlockSize, 0.0);
                    }
                }

                //reading dem buffer
                double xminBlock = xMin + xBlock*actualXBlockSize*pixelSpacing;
                double ymaxBlock = yMax - yBlock*actualYBlockSize*pixelSpacing;

                double xmaxBlock = xMin + (xBlock+1)*actualXBlockSize*pixelSpacing;
                double yminBlock = yMax - (yBlock+1)*actualYBlockSize*pixelSpacing;

                int cols = round((xmaxBlock - xminBlock)/demGeoTransform[1]) + 3;
                int rows = round((ymaxBlock - yminBlock)/demGeoTransform[1]) + 3;

                float *inDEMDataBuffer;
                inDEMDataBuffer = new float[cols*rows];

                int colOffset = int(-(demGeoTransform[2]*(demGeoTransform[3] - ymaxBlock) - demGeoTransform[0]*demGeoTransform[5] + demGeoTransform[5]*xminBlock)/(demGeoTransform[2]*demGeoTransform[4] - demGeoTransform[1]*demGeoTransform[5]));
                int rowOffset = int( (demGeoTransform[1]*(demGeoTransform[3] - ymaxBlock) - demGeoTransform[0]*demGeoTransform[4] + demGeoTransform[4]*xminBlock)/(demGeoTransform[2]*demGeoTransform[4] - demGeoTransform[1]*demGeoTransform[5]));

                inDEMBand->AdviseRead(colOffset, rowOffset, cols, rows, cols, rows, GDT_Float32, NULL);
                inDEMBand->RasterIO(GF_Read, colOffset, rowOffset, cols, rows, inDEMDataBuffer, cols, rows, GDT_Float32, 0, 0);

                for (int i = 0; i < actualXBlockSize; i++) { //for each col in block
                    x = xMin + (xBlock*actualXBlockSize + i + 0.5)*pixelSpacing;
                    for (int j = 0; j < actualYBlockSize; j++) {//for each row in block
                        y = yMax - ( actualYBlockSize*yBlock + j + 0.5)*pixelSpacing;
                        double height = getHeight(x, y, cols, colOffset, rowOffset, inDEMDataBuffer);

                        if (height > 0.0) {
                            OGRPoint *earthPoint;
                            //cout <<x <<"\t" <<y << endl;
                            earthPoint = new OGRPoint(x, y, height);
                            earthPoint->transform(transfTargetGeocentric);
                            OrbitData *sensorPositionFeature;
                            sensorPositionFeature = new OrbitData();
                            x = earthPoint->getX();
                            y = earthPoint->getY();
                            height = earthPoint->getZ();
                            if ( getPosition(x, y, height, sensorPositionFeature) )
                                for (auto& band : inputBands )
                                    bandBuffer[ band.first ][actualXBlockSize*j + i ] = getPixelValue(sensorPositionFeature, band.first);

                            delete sensorPositionFeature;
                            OGRGeometryFactory::destroyGeometry(earthPoint);
                        }
                    }
                }

                delete[] inDEMDataBuffer;
                inDEMDataBuffer = nullptr;

                for (auto& band : inputBands )
                    outBands[band.first]->WriteBlock(xBlock, yBlock, bandBuffer[band.first]);



            }
        }

        //clearing memory
        GDALClose(inDEMData2);
        GDALClose(outputDataset);
        /*
        for (auto& dataset:outDatasets)
            GDALClose(dataset.second);
        */
    }
}

double RangeDopplerGeocoding::getEarthPointZeroDopplerTime(double &x, double &y, double &z) {
    int lowerBound = 0;
    int upperBound = orbit->sensorPositionVectors.size() - 1;

    double lowerBoundFreq = getDopplerFrequency(x, y, z, lowerBound);
    double upperBoundFreq = getDopplerFrequency(x, y, z, upperBound);

    if (lowerBoundFreq == 0.0)
        return firstLineTime + lowerBound * lineTimeInterval;
    else if (upperBoundFreq == 0.0)
        return firstLineTime + upperBound * lineTimeInterval;
    else if (lowerBoundFreq * upperBoundFreq > 0.0)
        return nullDopplerTime;

    //starting binary search

    double midFreq = 0.0;

    //binary search to find the zero doppler time

    int mid = 0;
    while (upperBound - lowerBound > 1) {
        mid = (int) (lowerBound + upperBound) /2.0;

        OrbitData tmpFeature = orbit->sensorPositionVectors[mid];
        midFreq = tmpFeature.x_vel*(x - tmpFeature.x_pos) +
                tmpFeature.y_vel*(y - tmpFeature.y_pos) +
                tmpFeature.z_vel*(z - tmpFeature.z_pos);

        if (midFreq * lowerBoundFreq > 0.0) {
            lowerBound = mid;
            lowerBoundFreq = midFreq;
        }
        else if (midFreq * upperBoundFreq > 0.0) {
            upperBound = mid;
            upperBoundFreq = midFreq;
        }
        else if (midFreq == 0.0)
            return firstLineTime + mid * lineTimeInterval;
        //cout << upperBound << "\t" << lowerBound <<"\t" << midFreq << endl;
    }
    double y0 = lowerBound - lowerBoundFreq * (upperBound - lowerBound) / (upperBoundFreq - lowerBoundFreq);
    return firstLineTime + y0 * lineTimeInterval;
}

bool RangeDopplerGeocoding::getPosition(double &x, double &y, double &z, OrbitData *sensorPositionFeature) {
    double zeroDopplerTime = getEarthPointZeroDopplerTime(x, y, z);

    if (zeroDopplerTime == nullDopplerTime)
        return false;


    sensorPositionFeature->slant_range = computeSlantRange(zeroDopplerTime, x, y, z, sensorPositionFeature);

    if (!skipBistaticCorrection) { // skip bistatic correction for COSMO, TerraSAR-X and RadarSAT-2
        zeroDopplerTime += sensorPositionFeature->slant_range / lightSpeedInMetersPerDay;
        sensorPositionFeature->slant_range = computeSlantRange(zeroDopplerTime, x, y, z, sensorPositionFeature);
    }

    sensorPositionFeature->range_index = computeRangeIndex(zeroDopplerTime, sensorPositionFeature->slant_range);

    if (sensorPositionFeature->range_index < 0.0)
        return false;

    if (!nearRangeOnLeft)
        sensorPositionFeature->range_index = sourceImageWidth - 1 - sensorPositionFeature->range_index;

    sensorPositionFeature->azimuth_index = (zeroDopplerTime - firstLineTime)  / lineTimeInterval;

    //checking if cell is valid - perhaps it is required to be extended
    if ( sensorPositionFeature->range_index < 0 || sensorPositionFeature->range_index >= sourceImageWidth - 1 || sensorPositionFeature->azimuth_index > sourceImageHeight - 1 )
        return false;
    return true;
}

RangeDopplerGeocoding::RangeDopplerGeocoding(char *inFile, char *inDEM, int targetsrsCode, char *outFile): fileName(inFile), inDEMFileName(inDEM),
    skipBistaticCorrection(false), srgrFlag(true), isPolsar(true), nearEdgeSlantRange(0.0), sourceImageHeight(0), xBlockSize(0), yBlockSize(0),
    sourceImageWidth(0), sourceSRSEPSG(4326), targetSRSEPSG(targetsrsCode), nearRangeOnLeft(true), outputFile(outFile) {
    GDALAllRegister();

    sourceSRS.importFromEPSG(sourceSRSEPSG);
    targetSRS.importFromEPSG(targetSRSEPSG);

    metaDataFile = xmlParseFile(inFile);
    inPath = fileName.substr(0, fileName.length()-3)+"data/";

    //retrieving metadata
    parseMetaDataFile();

    //reading datasets
    readBands();

    //creating orbit info
    orbit = new Orbit(orbitVector, firstLineTime, lineTimeInterval, sourceImageHeight);
}

RangeDopplerGeocoding::~RangeDopplerGeocoding() {
    delete orbit;
    orbit = nullptr;

    delete demGeoTransform;
    demGeoTransform = nullptr;
}

void RangeDopplerGeocoding::setTargetExtents(double &xmin, double &ymin, double &xmax, double &ymax) {
    xMin = xmin;
    yMin = ymin;
    xMax = xmax;
    yMax = ymax;
}
