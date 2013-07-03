/* Interface for each lens element
   Computes the refraction at every surface
*/

#ifndef LENS_H
#define LENS_H

#include <iostream>
#include <fstream>
#include "geometry.h"
#include "loadData.h"

class Lens{
public:
	full_lens_data *data;
	ofstream out_file;
	Lens(full_lens_data *lens_data); // constructor
	float ComputeRefractedRay(Ray* incident_ray , Ray* refracted_ray, int lens_number, Point *Porigin);//, ofstream outFile);
private:
	int RaySphere(Point p1,Point p2,Point sc,double r,double *mu1,double *mu2); 
	int isInsideAperture(Point intersection, float aperture);
};

#endif
