/* Interface for each lens element
   Computes the refraction at every surface
*/

#include "lens.h"

#include <math.h>

using namespace std;

Lens::Lens(full_lens_data *lens_data){
	data = new full_lens_data;
	data = lens_data;
}

float Lens::ComputeRefractedRay(Ray* incident_ray, Ray* refracted_ray, int lens_number, Point *Porigin){
	lens_el *lensData = new lens_el;
	lensData = (*data)[lens_number];

	float index1 = lensData->ref_index;
	float index2;
	if(lens_number > 0){
		index2 = (*data)[lens_number-1]->ref_index;
	}else{
		index2 = 1;
	}
	double radius = (double) lensData->radius;
	float z_pos = lensData->z_pos;
	float aperture = lensData->aperture;

	Point sphere_center = Point(0,0,z_pos - radius);

	Vector normal = Vector(0,0,-1);

	float t = -1*(-1*incident_ray->o.z + (*data)[lens_number]->z_pos)/(Dot(incident_ray->d, normal));

	Point newIntPoint = incident_ray->o + t*incident_ray->d;

	//cout << newIntPoint.z << " " << "newPoint" << endl;

	Point Pcamera = incident_ray->o;
	
	double mu1;
	double mu2;

	Point intersectionPoint;

	if(radius != 0.){ // To avoid computing intersection of sphere and ray when in the center (aperture stop)
		int intersection = RaySphere(Pcamera,newIntPoint,sphere_center,fabs(radius),&mu1,&mu2);

		Point int_point1;
		Point int_point2;

		if(intersection){
			int_point1 = Pcamera + (mu1)*(newIntPoint - Pcamera);
			int_point2 = Pcamera + (mu2)*(newIntPoint - Pcamera);
		
		}else{
			//cerr <<"intersection could not be computed" << endl;
			return 0.f;
		}
	
		if(fabs(int_point1.z - z_pos) > fabs(int_point2.z - z_pos)){
			intersectionPoint = int_point2;
		}else{
			intersectionPoint = int_point1;
		}		
		
		if(!(fabs(intersectionPoint.x) < aperture/2 && fabs(intersectionPoint.y) < aperture/2)){
			return 0.f; // if ray hits outside the aperture return zero weight
		}
	
		//Ray *ray = new Ray();
		//ray->o = Pcamera;
		//ray->d = Normalize(intersectionPoint - Pcamera);
		Ray *ray = incident_ray;
		Vector normal_vector;
		//if(radius < 0.){
	//		normal_vector = Normalize(intersectionPoint - sphere_center);
	//	}else{
			normal_vector = Normalize(sphere_center - intersectionPoint);
	//	}

		//Refracted ray	
		Vector in_ray = -1*ray->d;
		float n = index1/index2;
		float c1 = -1*Dot(in_ray, normal_vector);
		float c2 = 1 - n*n*(1-c1*c1);
		Vector ref_ray;
		if(c2 > 0.0){
			c2 = sqrt(c2);
			ref_ray = Normalize(-1*(n*in_ray + (n*c1 - c2)*normal_vector));
			(*refracted_ray).o = intersectionPoint;
			(*refracted_ray).d = ref_ray;
		}else{
			//cerr << "Cannot compute refracted ray" << n << endl;
			return 0.f; 
		}
	}else{
		intersectionPoint = newIntPoint;
		*Porigin = newIntPoint;
		if(!(fabs(intersectionPoint.x) < aperture/2 && fabs(intersectionPoint.y) < aperture/2)){
			return 0.f; // if ray hits outside the aperture return zero weight
		}
		(*refracted_ray).o = intersectionPoint;//incident_ray->o;
		(*refracted_ray).d = incident_ray->d;
	}
	//outFile << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z;
	
	//cout << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << endl;
	return 1.f;	
}

int Lens::RaySphere(Point p1,Point p2,Point sc,double r,double *mu1,double *mu2){
   double a,b,c;
   double bb4ac;
   Point dp = Point(0,0,0);

   dp.x = p2.x - p1.x;
   dp.y = p2.y - p1.y;
   dp.z = p2.z - p1.z;
  
   a = dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
   b = 2 * (dp.x * (p1.x - sc.x) + dp.y * (p1.y - sc.y) + dp.z * (p1.z - sc.z));
   c = sc.x * sc.x + sc.y * sc.y + sc.z * sc.z;
   c += p1.x * p1.x + p1.y * p1.y + p1.z * p1.z;
   c -= 2 * (sc.x * p1.x + sc.y * p1.y + sc.z * p1.z);
   c -= r * r;
   bb4ac = b * b - 4 * a * c;
   if (fabs(a) == 0 || bb4ac < 0) {
      *mu1 = 0;
      *mu2 = 0;
      return 0;
   }
   
   *mu1 = (-b + sqrt(bb4ac)) / (2 * a);
   *mu2 = (-b - sqrt(bb4ac)) / (2 * a);
   
   return 1;
}

int isInsideAperture(Point intersection, float aperture){
	if(fabs(intersection.x) < aperture/2 && fabs(intersection.y) < aperture/2){
		return 1;
	}else{
		return 0;
	} 
}

int RaySphereIntersection(Ray *ray, Point sc, double r, double *mu1,double *mu2){

}
