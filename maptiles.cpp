/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */

      vector<Point<3>> points;
      map<Point<3>, TileImage*> createdMap;

      typename std::vector<TileImage>::iterator it;
      for(it = theTiles.begin(); it != theTiles.end(); it++) {
        LUVAPixel current = it->getAverageColor();
        Point<3> curr_point = convertToXYZ(current);
        points.push_back(curr_point);
        createdMap[curr_point] = &*it;
      }

      KDTree<3>* canvas = new KDTree<3>(points);
      MosaicCanvas* myCanvas = new MosaicCanvas(theSource.getRows(), theSource.getColumns());

      for(int i = 0; i < theSource.getRows(); i++) {
        for(int j = 0; j < theSource.getColumns(); j++) {
          LUVAPixel region = theSource.getRegionColor(i, j);
          Point<3> regionColor = convertToXYZ(region);
          Point<3> compare = canvas -> findNearestNeighbor(regionColor);
          TileImage* closest = createdMap[compare];
          myCanvas -> setTile(i, j, closest);
        }
      }

      delete canvas;
      canvas = NULL;
      return myCanvas;
}
