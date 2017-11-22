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

int main(int argc, char *argv[]) {

    if (argc < 6) {
        cout <<"usage: ./RangeDopplerTerrainCorrection sentinel_metatada.dim dem_file EPSG pixelSize output_directory xMin yMin xMax yMax\n"
               "xMin yMin xMax yMax: optional\n";
        return -1;
    }
    cout << "Starting Range-Doppler Geocoding\n";
    //cout << setprecision(20) << endl;
    RangeDopplerGeocoding geocode(argv[1], argv[2], atoi(argv[3]),  argv[5]);
    //get tie point grid and compute new image dimensions
    if (argc == 6)
        geocode.getTiePointGrid();
    else if (argc == 10) {
        double xMin = atof(argv[6]);
        double yMin = atof(argv[7]);
        double xMax = atof(argv[8]);
        double yMax = atof(argv[9]);

        geocode.setTargetExtents(xMin, yMin, xMax, yMax);

    }
    //resample image
    float pixelSize = atof(argv[4]);
    geocode.resampleImage(pixelSize);

    cout << "All done!\n";

    return 0;
}
