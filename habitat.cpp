/****************************  habitat.cpp   **********************************
* Author:        Agner Fog
* Date created:  2023-09-10
* Last modified: 2024-10-13
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines functions for manipulating territories
* and moving boudaries between territories on a square piece of land
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"


void SortedPointList::put(TPoint const & p) {    // add a point to list
    // binary search in list
    int a = 0, b = npoints, i;
    while (b > a) {
        i = (a + b) / 2;
        if (p < pointlist[i]) b = i; else a = i + 1;
    }
    i = a - 1;
    // check if p already in list
    if (i >= 0 && pointlist[i] == p) return;
    // check if list is full
    if (npoints >= pointListSize) {
        errors.reportError("Error in function SortedPointList::put: List overflow");
        return;
    }
    // put point in list at position i+1
    if (npoints > a) {
        memmove(pointlist + a + 1, pointlist + a, (npoints - a) * sizeof(TPoint));
    }
    pointlist[a] = p; npoints++;
}

TPoint SortedPointList::pop(void) {    // get lowest point out of list
    if (npoints == 0) return TPoint(0, 0);
    TPoint p = pointlist[0];
    npoints--;
    memmove(pointlist, pointlist + 1, npoints * sizeof(TPoint));
    return p;
}


void SortedPDList::put(TPoint p) {
    TPointDist q = makePoint(p);
    // binary search in list
    int a = 0, b = npoints, i;
    while (b > a) {
        i = (a + b) / 2;
        if (q < pointlist[i]) b = i; else a = i + 1;
    }
    i = a - 1;
    // check if q already in list
    if (i >= 0 && pointlist[i] == q) return;
    // check if list is full
    if (npoints >= pointListSize) {
        errors.reportError("Error in function SortedPDList::put: List overflow");
        return;
    }
    // put point in list at position i+1
    if (npoints > a) {
        memmove(pointlist + a + 1, pointlist + a, (npoints - a) * sizeof(TPointDist));
    }
    pointlist[a] = q;
    npoints++;
}

TPoint SortedPDList::pop(void) {
    // get point with highest ddist from list
    if (npoints == 0) return TPoint(0, 0);
    else return pointlist[--npoints];  // convert to base class TPoint (ignore "do not slice" warning)
}

void SortedPDList::remove(TPoint p) {  // remove a specific point from list
    if (npoints == 0) return;
    int i, b;
    i = find(p);
    if (i < 0) return;
    b = i + 1;
    if (npoints - b > 0) {
        memmove(pointlist + i, pointlist + b, (npoints - b) * sizeof(TPointDist));
    }
    npoints--;
}

int SortedPDList::find(TPoint p) {
    // find a point in list and return its index
    // (returns -1 if not found)
    TPointDist q = makePoint(p);
    // binary search in list
    int a = 0, b = npoints, i;
    while (b > a) {
        i = (a + b) / 2;
        if (q < pointlist[i]) b = i; else a = i + 1;
    }
    i = a - 1;
    if (i >= 0 && pointlist[i] == q) {
        return i;        // point found
    }
    return -1;           // point not found
}

// member function for class Habitat

Habitat::Habitat(AltruData * d) {                             // constructor
    setData(d);
}

void Habitat::setData(AltruData * d) {                        // set all data
    // set parameters for area
    this->d = d;
    map = (idType*)(d->extraBuffer[ebMap]);                   // memory allocated for map of points
    buffer1 = (TPoint*)d->extraBuffer[ebPointsSorted];        // memory allocated for lists of points
    buffer2 = (TPointDist*)d->extraBuffer[ebPointDistSorted]; // memory allocated for lists of points    
    rowLength = d->rowLengthTerri;
    numRows = d->numRows;
}

// Find an arbitrary point on the border of a territory.
// The return value is a point on the border
TPoint Habitat::findBorder(TerriGroup const * d1) const {
    if (d1->id == 0) return TPoint(0, 0);
    // find a point that is mine
    TPoint p0(d1->centerx, d1->centery);
    TPoint p1;
    idType d2 = owner(p0);
    int dir = 0;
    if (d2 != d1->id) {
        // territory does not own its own center of gravity
        p0 = findCenter(d1);
    }
    // find a border point
    p1 = p0;
    do {
        p0 = p1;
        p1 = step(p0, dir);
    } while (owner(p1) == d1->id);
    return p0;
}


// Find an arbitrary point on the border of a territory.
// The return value is a point on the border
// This is the same as findBorder, but this version tries to avoid border to an enclave.
// This is pretty slow, but necessary for graphics output
TPoint Habitat::findBorder2(TerriGroup const * d1) const {
    if (d1->id == 0) return TPoint(0, 0);
    // find a point that is mine
    int dir = 0;
    TPoint p0(d1->centerx, d1->centery);
    TPoint p1, p2;
    if (owner(p0) != d1->id) {
        // territory does not own its own center of gravity. may be an enclave
        p0 = findCenter(d1);
    }
    int32_t biggestPerimeter = 0;
    idType  neighb = 0;
    idType  firstNeighb = 0;
    idType  lastNeighb = 0;
    TPoint  bestBorderPoint = TPoint(0,0);

    // search in all directions for border points
    for (dir = 0; dir <= 3; dir++) {
        p1 = p0;
        while (true) {
            p2 = step(p1, dir);
            neighb = owner(p2);
            if (neighb != d1->id) break;
            p1 = p2;
        }
        if (neighb > 0 && neighb != firstNeighb && neighb != lastNeighb) {
            // found a new neighbor
            // measure perimeter of own territory, starting at border to this neighbor
            int perimeter = 0;
            borderWalk(d1->id, p1, [&perimeter](TPoint&){
                perimeter++;
                }, 
                [](TPoint&,TPoint&){}
            );
            if (perimeter > biggestPerimeter) {
                biggestPerimeter = perimeter;
                bestBorderPoint = p1;
            }
        }
        lastNeighb = neighb;        
        if (dir == 0) firstNeighb = neighb;
    }
    return bestBorderPoint;
}


// findNeighbors
// This function finds all neighbors to a territory and returns information
// needed by the function conquer.
// Information about each neighbor is returned in the list neighborList.
// Each record contains the id of the territory, the first border point, the length of
// the shared border, and information on whether the border is contiguous.
// The return value is the number of records in neighborList.

// The function follows the border of the territory d1 all way round in search
// for neighbors. It also traverses the territory in the east-west and
// north-south directions to search for any other territory that is completely
// surrounded by d1 (enclave). This method is not 100% safe: There remains a
// possibility that a small enclave surrounded by d1 remains undetected.
// In the extremely rare case of two weirdly shaped enclaves or one very
// convoluted enclave, this algorithm may find the enclaves only.
// Eliminating these possibilities of failure completely would require that
// all points in the entire area be checked, which would be extremely
// time-consuming. The possible failure will be extremely rare and without
// serious consequences.
int Habitat::findNeighbors(TerriGroup const * const d1, Neighbor * const neighborList, int * const numNeighbors) const {
    int i;                        // loop counter
    TPoint p0;                    // starting point
    TPoint p3;                    // search for border, start of section
    TPoint p4;                    // center of territory d1
    idType d3;                    // last neighbor id
    if (d1->id == 0) return 0;
    *numNeighbors = 0;            // reset number of neighbors
    if (d1->area == 0) return 0;

    p4 = TPoint(d1->centerx, d1->centery);       // center of territory
    idType d2 = owner(p4);
    int dir1 = 0;                                // direction
    if (d2 != d1->id) {
        // territory does not own its own center of gravity (this is rare)
        p4 = findCenter(d1);
    }
    // search in all four directions from p4 for border and enclaves
    for (dir1 = 0; dir1 < 4; dir1++) {
        p3 = p4;
        do {
            p0 = p3;
            p3 = step(p3, dir1);
            d2 = owner(p3);
        } while (d2 == d1->id);
        // now p0 is a border point. Check if this neighbor is already in list
        for (i = 0; i < *numNeighbors; i++) {
            if (neighborList[i].id == d2) {
                i = -1; break;
            }
        }
        if (i == -1) continue;                   // already found this border. continue with next dir1

        // walk around newly found border to find all border points
        int nn1 = *numNeighbors;                 // number of neighbors before border walk
        d3 = d1->id;
        borderWalk(d1->id, p0, [](TPoint&) {}, [this, numNeighbors, neighborList, &d3, &p3](TPoint&p_1, TPoint&p_2) {
            // do_them lambda function:
            // check if neighbor is in list
            int i = 0;
            int32_t owner2 = owner(p_2);
            if (owner2 > 0) {
                for (i = *numNeighbors - 1; i >= 0; i--) {
                    if (neighborList[i].id == owner2) {
                        // update record
                        neighborList[i].sharedBorder++;
                        if (i < *numNeighbors - 1) {
                            // not contiguous
                            if (owner2 != d3) {
                                neighborList[i].nonContiguous++;
                                d3 = owner2;
                                //??p3 = p_1;
                                p3 = p_2;
                            }
                        }
                        break;
                    }
                }
            }
            if (i < 0) {
                // this neighbor is not in list. add to list
                if (*numNeighbors < maxNeighbors) {
                    i = (*numNeighbors)++;
                    neighborList[i].id = owner(p_2);
                    //??neighborList[i].start = p_1;
                    neighborList[i].start = p_2;
                    neighborList[i].sharedBorder = 1;
                    neighborList[i].nonContiguous = 0;
                }
                else {
                    errors.reportError("Error in function Habitat::findNeighbors: Neighbor list overflow");
                    (*numNeighbors)--;
                }
            }
        }); // end of BorderWalk 
        if (neighborList[nn1].nonContiguous && d3 == neighborList[nn1].id) {
            // join start and end point into contiguous section
            neighborList[nn1].nonContiguous--;
            neighborList[nn1].start = p3;
        }
    }
    return *numNeighbors;
}

// check of d1 and d2 are neighbors
// returns a point on neighbor border if true, TPoint(0,0) if not
TPoint Habitat::isNeighbor(TerriGroup const * d1, idType d2) const {
    bool found = false;
    TPoint p1 = findBorder(d1);
    TPoint pn(0, 0);  // neighbor point
    borderWalk(d1->id, p1, [](TPoint&) {}, [this, d2, &found, &pn](TPoint & p_1, TPoint & p_2) {
        // do_them lambda function
        if (!found && owner(p_2) == d2) {
            found = true;
            pn = p_2;
        }        
    });
    return pn;
}


// transfer na units of land from owner 'nfrom' to conqueror 'to'
void Habitat::conquer(TerriGroup * from, TerriGroup * to, int32_t na, Neighbor const * nfrom) const {
    // Transfer 'na' bits of land from territory 'from'
    // to 'to'. 'nfrom' is information returned by
    // function findNeighbors on the 'from' territory.

    if (from->id == to->id) return;

    // make list of points on 'from's border to 'to'
    TPoint centerFrom(from->centerx, from->centery);
    TPoint centerTo(to->centerx, to->centery);

    SortedPDList list(centerFrom, centerTo, buffer2);
    TPoint p0, p1, p2, p3;       // points on border
    int dir;                     // direction
    int breaks;                  // number of discontinuities in border to 'to'
    int np2 = 0;                 // count borderpoints
    int32_t tarea = 0;           // transferred area
    uint32_t sx = 0, sy = 0;     // sum of x,y coordinates of transferred area
    idType   newid = to->id;

    if (na >= from->area - d->territorySizeMin) {
        // take it all
        TPoint tb = findBorder(from);
        floodFill(tb, to->id);
        return;    
    }

    if (nfrom) {
        // normal call
        breaks = nfrom->nonContiguous;
        p0 = nfrom->start;

        if (owner(p0) != from->id) {
            errors.reportError("Error in function Habitat::conquer: point not owned");
            return;
        }
        bool stop = false;
        // walk around border of 'from' to find points bordering to NewOwner
        TPoint pb = findBorder(from);
        borderWalk(from->id, pb, [](TPoint&) {}, 
            [this, newid, &list, &np2, &breaks, &stop](TPoint& p_1, TPoint& p_2) {
                // do_them lambda function            
                if (!stop) {
                    if (owner(p_2) == newid) {
                        list.put(p_1);
                    }
                    else {
                        if (list.npoints > np2) {
                            np2 = list.npoints;
                            if (--breaks < 0) stop = true;
                        }
                    }
                }
            }
        );
    }
    else {
        // called from splitArea (centerTo is not yet owned by 'to')
        changeOwner(centerTo, newid);
        list.put(centerTo);
    }

    if (list.npoints == 0) {
        // 'from' is not a neighbor. may have been lost since findNeighbors was called
        return;
    }

    // transfer points
    while (tarea < na && list.npoints) {
        p1 = list.pop();  // get nearest point
        if (owner(p1) != from->id && nfrom) {
            errors.reportError("Error in function Habitat::conquer: point not owned");
        }
        changeOwner(p1, to->id);
        tarea++;
        sx += p1.x; sy += p1.y;
        // check neighbor points
        for (dir = 0; dir < 4; dir++) {
            p2 = step(p1, dir);
            if (owner(p2) == from->id) {
                list.put(p2);        // add new border point to list
            }
        }
    }
   

    // update information in groups after transfer
    updateTerr(from, to, tarea, sx, sy, d->carryingCapacity[0]);

    if (tarea < na) {
        // failed to transfer area
        errors.reportError("Error in function Habitat::conquer: Failed to transfer area");
    }

    if (list.npoints == 0) return;

    // Check if oldOwner territory is still contiguous by walking around
    // the border of oldOwner and see if we can find all the points in list:
    p0 = list.get(list.npoints - 1);
    p3 = TPoint(0, 0);
    bool stop = false;
    borderWalk(from->id, p0, [](TPoint&) {}, [this, to, &p3, &list, &stop](TPoint & p_1, TPoint & p_2) {
        // do_them lambda function
        if (!stop) {
            if (owner(p_2) == to->id && p_1 != p3) {
                p3 = p_1; // avoid trying same point twice
                list.remove(p3);
                if (list.npoints == 0) stop = true;
            }
        }
    });

    if (list.npoints > 0) {
        // Didn't find all points in list, oldOwner territory has become noncontiguous.
        // Find all sections of OldOwner territory
        list.put(p0); // put p0 back in list to represent the section we have found.

        // make list of which section each point belongs to
        // (this happens so rarely that we can afford not to recycle allocated memory)
        int32_t * pslist = new int32_t[list.npoints+1]();

        int i, k, L, nsect = 0, perimeter, nbig = 0, lbig = 0;
        for (i = 0; i < list.npoints; i++) {
            if (pslist[i] == 0) {
                // found a section
                pslist[i] = ++nsect;
                // walk around boundary of section to measure its perimeter
                p0 = list.get(i);
                perimeter = 0;
                borderWalk(from->id, p0, [&list, &pslist, nsect, &perimeter](TPoint & p_1) {
                    int j = list.find(p_1);
                    if (j >= 0) 
                        pslist[j] = nsect;
                    perimeter++;
                    }
                , [](TPoint&,TPoint&) {});
                if (perimeter > lbig) {
                    // biggest section so far
                    lbig = perimeter; nbig = nsect;
                }
            }
        }

        // Now nsect is the number of sections the OldOwner territory is divided into.
        // nbig identifies the biggest section, lbig is its length.
        // Let NewOwner steal all but the biggest section
        tarea = 0;
        for (k = 1; k <= nsect; k++) {
            if (k != nbig) {
                // find one point belonging to this section
                for (L = 0; L < list.npoints; L++) {
                    if (pslist[L] == k) {
                        TPoint pl = list.get(L);
                        if (owner(pl) == from->id) { // check if point still belongs to old owner //??
                            tarea += floodFill(pl, to->id);
                        }
                        break;
                    }
                }
            }
        }        
        delete[] pslist;
    }
}


// change owner of all points of land contiguous with p        
int32_t Habitat::floodFill(TPoint const & p, idType newOwner) const {

    // change owner of all points contiguous with p
    if (owner(p) == newOwner) {
        errors.reportError("Error in function Habitat::floodFill: same owner");
        return 0;
    }
    SortedPointList list(buffer1);               // sorted list of points to change owner of
    TPoint p0;                                   // point from list
    TPoint p1;                                   // point to put in list
    int dir;
    int32_t oldOwner = owner(p);
    uint32_t tarea = 0;                          // transferred area
    uint32_t sx = 0, sy = 0;                     // sum of x,y coordinates of transferred area
    list.put(p);                                 // put first point in list
    while (list.npoints) {                       // continue as long as there are points in list
        p0 = list.pop();                         // get point from list
        changeOwner(p0, newOwner);               // change owner of point
        tarea++;                                 // update statistics
        sx += p0.x;
        sy += p0.y;
        // check all neighbor points to p0
        for (dir = 0; dir < 4; dir++) {
            p1 = step(p0, dir);
            if (owner(p1) == oldOwner) {
                // put p1 into list if it is not already there
                list.put(p1);
            }
        }
    }

    // update information in groups after transfer
    TerriGroup * from = (TerriGroup*)(d->groupData) + oldOwner;
    TerriGroup * to   = (TerriGroup*)(d->groupData) + newOwner;
    updateTerr(from, to, tarea, sx, sy, d->carryingCapacity[0]);
    return (int32_t)tarea;
}


// split territory into two, the second with area ar
void Habitat::splitArea(TerriGroup * oldOwner, TerriGroup * newOwner, int32_t ar) const {

    // split territory into two. area of new territory is ar
    if (ar >= oldOwner->area || int32_t(ar) <= 0) return;
    if (oldOwner->area < 2 * ar) ar = oldOwner->area - ar;
    if (newOwner->area) {
        errors.reportError("Error in function Habitat::splitArea: NewOwner not empty");
        return;
    }

    uint32_t sdist = 0;
    TPoint p0 = findBorder(oldOwner);
    TPoint p3 = p0;
    TPoint center(oldOwner->centerx, oldOwner->centery); // center of old owner
    sdist = 0;

    // walk around boundary to find border point p3 with highest distance from center
    borderWalk(oldOwner->id, p0, [center, &sdist, &p3](TPoint & p_1) {
        // do_us lambda function
        uint32_t dd = squareDistance(p_1, center);
        if (dd > sdist) {
            sdist = dd; p3 = p_1;
        }
        },
        [](TPoint&, TPoint&) {}
    );

    // point p3 is the most distant. Make it a base for new territory
    newOwner->centerx = p3.x;  newOwner->centery = p3.y;  
    newOwner->sx = newOwner->sy = 0; newOwner->area = 0; newOwner->age = 0;

    // grow the new territory
    conquer(oldOwner, newOwner, ar, 0);
}


// update information in both territories after transfer of land
void Habitat::updateTerr(TerriGroup * oldOwner, TerriGroup * newOwner, int32_t tarea, uint32_t sx, uint32_t sy, float popDens) const {
    // update information in both territories after transfer
    oldOwner->sx -= sx; oldOwner->sy -= sy;
    oldOwner->area -= tarea;
    if (oldOwner->area > 0) {
        oldOwner->centerx = oldOwner->sx / oldOwner->area;
        oldOwner->centery = oldOwner->sy / oldOwner->area;
    }
    else if (oldOwner->area == 0) {
        //newOwner->nn += oldOwner->nn;
        oldOwner->nn = 0;
        oldOwner->centerx = 0;  oldOwner->centery = 0;
    }
    else {
        errors.reportError("negative area after conquer");
    }
    newOwner->sx += sx; newOwner->sy += sy;
    newOwner->area += tarea;
    if (newOwner->area) {
        newOwner->centerx = newOwner->sx / newOwner->area;
        newOwner->centery = newOwner->sy / newOwner->area;
    }
    else {
        newOwner->nn = 0;
        newOwner->centerx = 0;  newOwner->centery = 0;
    }
    oldOwner->nmax = (int32_t)lround(oldOwner->area * popDens);
    newOwner->nmax = (int32_t)lround(newOwner->area * popDens);
}

// find a point belonging to a territory in the rare case that 
// it does not own its own center of gravity
TPoint Habitat::findCenter(TerriGroup const * d1) const {
    if (d1->area == 0) return TPoint(0,0);
    int32_t x, x1, x2, y, y1, y2;                // coordinates
    int32_t radius;
    idType  owner = d1->id;

    for (radius = 1; radius < rowLength; radius++) {
        // search in ever wider sqaures
        x1 = d1->centerx - radius; if (x1 < 1) x1 = 1;
        x2 = d1->centerx + radius; if (x2 > rowLength-2) x2 = rowLength-2;
        y1 = d1->centery - radius; if (y1 < 1) y1 = 1;
        y2 = d1->centery + radius; if (y2 > numRows-2) y2 = numRows-2;
        for (x = x2; x >= x1; x--) {
            y = y2;
            if (map[x+y*rowLength] == owner) return TPoint(x,y); // found a point belonging to d1
            y = y1;
            if (map[x+y*rowLength] == owner) return TPoint(x,y);
        }         
        for (y = y2-1; y > y1; y--) {
            x = x2;
            if (map[x+y*rowLength] == owner) return TPoint(x,y);
            x = x1;
            if (map[x+y*rowLength] == owner) return TPoint(x,y);
        }
    }
    errors.reportError("Cannot find any point belonging to territory");
    return TPoint(0,0);  // not found
}

// check if a territory is contiguous. This is used for debugging only
void Habitat::checkIfContiguous(TerriGroup const * d1) const {
    if (d1->area == 0) return;
    bool * checkmap = new bool[d->totArea+1]();
    TPoint p = findCenter(d1);
    SortedPointList list(buffer1);               // sorted list of points to change owner of
    TPoint p0;                                   // point from list
    TPoint p1;                                   // point to put in list
    int dir;
    int32_t tarea = 0;                           // observed area
    list.put(p);                                 // put first point in list
    while (list.npoints) {                       // continue as long as there are points in list
        p0 = list.pop();                         // get point from list
        tarea++;                                 // count points
        checkmap[p0.x+p0.y*rowLength] = true;    // remember that this point has been counted
        // check all neighbor points to p0
        for (dir = 0; dir < 4; dir++) {
            p1 = step(p0, dir);
            if (owner(p1) == d1->id && !checkmap[p1.x+p1.y*rowLength]) {
                // put p1 into list if it has not already been counted
                list.put(p1);
            }
        }
    }
    if (tarea != d1->area) {
        char text[1024];  int lim = 0;
        sprintf(text, "Territory 0x%X is not contiguous in generation %i\nLost points:\n", d1->id, (int)d->generations);
        for (int i = 0; i < d->totArea; i++) {
            if (map[i] == d1->id && !checkmap[i] && ++lim < 50) { // found a lost point
                sprintf(text+strlen(text), "(%X,%X) ", i % d->rowLengthTerri, i / d->rowLengthTerri);            
            }        
        }
        errors.reportError(text);    
    }
    delete[] checkmap;
}

    
// check all territories for consistency. This is used for debugging only
void checkAllTerritories (AltruData * d) {
    int32_t * areaSums = new int32_t[d->nIslands+1]();
    int32_t * xSums    = new int32_t[d->nIslands+1]();
    int32_t * ySums    = new int32_t[d->nIslands+1]();
    idType  * map = (idType*)(d->extraBuffer[ebMap]); 

    for (int a = 0; a < d->totArea; a++) {
        int owner = map[a];
        if (owner > d->nIslands) {
            errors.reportError("owner id out of range");
            return;
        }
        if (owner > 0) {
            areaSums[owner]++;
            int x = a % d->rowLengthTerri;
            int y = a / d->rowLengthTerri;
            xSums[owner] += x;
            ySums[owner] += y;
        }
    }

    int iGroup;
    TerriGroup * group;
    for (iGroup = 1; iGroup <= d->nIslands; iGroup++) {
        group = (TerriGroup*)(d->groupData) + iGroup; // point to current group
        if (group->area != areaSums[iGroup]) {
            errors.reportError("Territory area wrong");
        }
        if (group->sx != xSums[iGroup] || group->sy != ySums[iGroup]) {
            errors.reportError("Territory coordinate sums wrong");
        }
        if (group->area > d->territorySizeMax * 4) {
            //errors.reportError("Territory area too big");         
        }
        if (group->area > 0 && group->area < d->territorySizeMin / 4) {
            //errors.reportError("Territory area too small");         
        }
    }
    delete[] areaSums;
    delete[] xSums;
    delete[] ySums;
}


void checkArea(AltruData * d) {
    // check and calculate area of terrain
    if (d->totArea <= 0) {
        d->totArea = 40000;
    }
    if (d->totArea != d->rowLengthTerri * d->numRows) {
        // area has changed. recalculate size
        d->rowLengthTerri = int(sqrtf(d->totArea));        // round down square root of area
        if (d->rowLengthTerri < 32) d->rowLengthTerri = 32;
        d->numRows = d->totArea / d->rowLengthTerri;
        if (d->numRows < 32) d->numRows = 32;
        d->totArea = d->rowLengthTerri * d->numRows;       // make sure area has no partial row at the end
        d->lastArea = d->totArea;                          // remember area
        if (d->territorySizeMin < 10) d->territorySizeMin = 10;
        if (d->territorySizeMax > d->totArea / 8) d->territorySizeMax = d->totArea / 8;
        d->maxIslands = d->totArea / d->territorySizeMin;  //  max number of groups
        d->minGroupSize = d->territorySizeMin * d->carryingCapacity[0] * 0.5f;
        if (d->minGroupSize < 4) d->minGroupSize = 4;
        if (d->maxIslands >= maxId) {                      // check limit if id is 16 bit number
            errors.reportError("Too many groups. Maximum is 65000");
            d->maxIslands = 65000;
            d->totArea = d->maxIslands * d->territorySizeMin;// recalculate
            d->rowLengthTerri = int(sqrtf(d->totArea));
            d->numRows = d->totArea / d->rowLengthTerri;
            d->totArea = d->rowLengthTerri * d->numRows;
            d->lastArea = d->totArea;
            d->minGroupSize = d->territorySizeMin * d->carryingCapacity[0] * 0.5f;
        }
    }
    d->nMaxPerGroup = (int32_t)lround(d->territorySizeMax * d->carryingCapacity[locusAltruism]);
}
