/****************************  habitat.h   ************************************
* Author:        Agner Fog
* Date created:  2023-09-10
* Last modified: 2023-12-31
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines common data structures, classes, and functions
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#pragma once

// constant definitions
int const ebMap = 0;                             // extraBuffer for map
int const ebPointsSorted = 1;                    // extraBuffer for sorted list of points
int const ebPointDistSorted = 2;                 // extraBuffer for sorted list of points
int const ebOrder = 3;                           // extraBuffer for random order list of groups
int const pointListSize = 4096;                  // maximum number of points in sorted point list
int const maxNeighbors = 32;                     // maximum number of neighbors in list

// type for group id. may use 16 bits to save memory in the huge bit point map
typedef uint16_t idType;                         // 16 bits
int const maxId = 0xFFFE;                        // maximum value of group id


// Deme structure, describing each deme or territory
struct TerriDeme {
    int32_t  nn;                                 // number of individuals in group
    int32_t  nmax;                               // carrying capacity
    int32_t  gAltruism[2];                       // gene pool, [0]: egoism or neutral; [1]: altruism or regality
    int32_t  area;                               // area of territory
    int32_t  emigrationPotential;                // excess available for emigration
    idType   id;                                 // index into demeData, used as owner for points belong to this deme
    idType   age;                                // group age
    uint32_t sx;                                 // sum of x coordinates
    uint32_t sy;                                 // sum of y coordinates
    uint16_t centerx, centery;                   // center of gravity
    float groupfit;                              // group fitness
};


// Structure and class definitions for manipulation of point map

// struct TPoint defines the x,y-coordinates of one point on the geographic map.
// Each point can belong to a group territory, found by Habitat::owner()
struct TPoint {
    uint16_t x;                                  // x coordinate
    uint16_t y;                                  // y coordinate
    TPoint() {};                                 // default constructor
    TPoint(uint16_t x, uint16_t y) {             // constructor
        this->x = x;  this->y = y;
    }
    uint32_t pack() const {
        return uint32_t(y) << 16 | x;            // pack into 32-bit integer
    }
    bool operator == (TPoint const & r) const {  // equal
        return pack() == r.pack();
    }
    bool operator != (TPoint const & r) const {  // not equal
        return pack() != r.pack();
    }
    bool operator < (TPoint const & r) const {   // less than
        return pack() < r.pack();
    }
};

static inline TPoint step(TPoint & p, int direction) {
    // find neighbor point
    switch (direction) {
    case 0: // east
        return TPoint(p.x + 1, p.y);
    case 1: // north
        return TPoint(p.x, p.y + 1);
    case 2: // west
        return TPoint(p.x - 1, p.y);
    case 3: // south
        return TPoint(p.x, p.y - 1);
    }
    return p;
}

static inline int32_t squareDistance(TPoint p1, TPoint p2) {
    // square distance between two points
    int32_t dx = int32_t(p1.x) - int32_t(p2.x);
    int32_t dy = int32_t(p1.y) - int32_t(p2.y);
    return dx * dx + dy * dy;
}

// define structure containing point and distance for use in SortedPDList
struct TPointDist : public TPoint {              // record in list: point and square distance
    int32_t ddist;                               // difference in square distance from two points
    TPointDist(void) {};                         // default constructor
    TPointDist(TPoint const & p) : TPoint(p) {}; //  constructor
    bool operator < (TPointDist const&r) const { // sort by distance
        return ddist < r.ddist || (ddist == r.ddist && pack() < r.pack());
    }
};

// define sorted list of TPoint used in function Habitat::splitArea

class SortedPointList {                          // sorted list of points for FloodFill
    TPoint * pointlist;                          // sorted list of points
public:
    int npoints;                                 // number of points in list
    SortedPointList(TPoint * buffer) {           // constructor
        npoints = 0;
        pointlist = buffer;
    }
    void put(TPoint const & p);                  // put a point in list
    TPoint pop(void);                            // get lowest point from list
};

// define sorted list of TPointDist used in function conquer
class SortedPDList {
public:
    SortedPDList(TPoint pfrom, TPoint pto, TPointDist * buffer) { // constructor
        // constructor
        npoints = 0; from = pfrom; to = pto;
        pointlist = buffer;
    }
    int npoints;                                 // number of points in list
    TPointDist makePoint(TPoint const & p) {     // convert TPoint to TPointDist
        TPointDist q(p);
        q.ddist = squareDistance(p, from) - squareDistance(p, to);
        return q;
    }
    void put(TPoint p);                          // put a point in list
    TPoint pop(void);                            // get point with highest ddist from list
    TPoint get(int i) {                          // get point with index i from list
        return pointlist[i];                     // convert to base class TPoint (ignore "do not slice" warning)
    }
    TPointDist getd(int i) {                     // get TPointDist with index i from list
        return pointlist[i];
    }
    void remove(TPoint p);                       // remove a specific point from list
    int find(TPoint p);                          // find point in list and return its index
protected:
    TPointDist * pointlist;                      // list of points
    TPoint from, to;                             // reference points for computing distance
};



// structure for list of neighbors to a given territory
struct Neighbor {                                // used for results from function findNeighbors
    int32_t sharedBorder;                        // length of shared border
    int32_t nonContiguous;                       // if border to this neighbor not contiguous
    TPoint  start;                               // first border point
    idType  id;                                  // neighbor deme
};


// class Habitat encapsulates manipulations on the geographic map.
// These functions are encapsulated in a class so that they have easy access to map, rowLength, etc.
class Habitat {
public:
    Habitat(AltruData * d);                      // constructor
    void setData(AltruData * d);                 // set all data
    idType owner(TPoint const &p) const {
        return map[p.x+p.y*rowLength];           // find owner of a point
    }
    void changeOwner(TPoint &p, idType newOwner) const { // change owner of a point
        if (p.x == 0 || p.y == 0) {
            errors.reportError("owning point beyond border");
        }

        map[p.x+p.y*rowLength] = newOwner;
    }
    // find a point on the border of a group territory
    TPoint findBorder(TerriDeme const * d1) const;
    // same as findBorder, tries to avoid border to an enclave:
    TPoint findBorder2(TerriDeme const * d1) const;
    // find all neighbor territories with a shared border with territory d1
    int findNeighbors(TerriDeme const * const d1, Neighbor * const neighborList, int * const numNeighbors) const;
    // transfer na units of land from owner 'nfrom' to conqueror 'to'
    void conquer(TerriDeme * from, TerriDeme * to, int32_t na, Neighbor const * nfrom) const;
    // change owner of all points of land contiguous with p        
    void splitArea(TerriDeme * oldOwner, TerriDeme * newOwner, int32_t ar) const;
    // check if two groups are neighbors
    TPoint isNeighbor(TerriDeme const * d1, idType d2) const;
    // find central point belonging to a territory
    TPoint findCenter(TerriDeme const * d1) const;
    // check if territory is contiguous. This is used for debugging only
    void checkIfContiguous(TerriDeme const * d1) const;
    
protected:
    int32_t floodFill(TPoint const & p, idType newOwner) const;
    // update information in both territories after transfer of land
    void updateTerr(TerriDeme * oldOwner, TerriDeme * newOwner, int32_t tarea, uint32_t sx, uint32_t sy, float popDens) const;
    AltruData * d;                               // pointer to common data structure
    idType * map;                                // pointer to map of point owners
    TPoint * buffer1;                            // memory allocated for list of points
    TPointDist * buffer2;                        // memory allocated for list of point distances
    int32_t   rowLength;                         // length of rows, size of terrain in the x direction
    int32_t   numRows;                           // number of rows, size of terrain in the y direction
public:
    // find all points on the borders between a territory and all its neighbors
    void borderWalk(idType us, TPoint & start, std::function<void(TPoint& p_1)> do_us, std::function<void(TPoint&p_1, TPoint&p_2)> do_them) const {
        // borderWalk is used by several functions to follow the border all the way around a territory.
        // 'us' is the owner of the territory, TPoint start is any point on the border,
        // do_us is a function to do for each border point p_1,
        // do_them is a function to do for each point p_2 on the foreign side of the border.
        // (this function is inlined for the sake of more efficient calling of do_us and do_them)
        if (owner(start) != us) {
            return;                              // error. should not happen
        }
        TPoint p_1 = start;                      // point on our border
        TPoint p_2;                              // point on their border
        int dir, dir0;                           // direction
        for (dir0 = 3; dir0 > 0; dir0--) {
            p_2 = step(p_1, dir0);
            if (owner(p_2) != us) break;
        }
        dir = dir0;
        do {
            p_2 = step(p_1, dir);
            if (owner(p_2) == us) {
                p_1 = p_2;
                dir = (dir + 1) & 3;
                do_us(p_1);
            }
            else {
                dir = (dir - 1) & 3;
                do_them(p_1, p_2);
            }
        } while (p_1 != start || dir != dir0);
    }
};

// check all territories for consistency. This is used for debugging only
void checkAllTerritories (AltruData * d);

// check and calculate area of terrain
void checkArea(AltruData * d);
