#include <iostream>
#include <map>

#ifndef LOADDATA_H
#define LOADDATA_H

using namespace std;

struct lens_element{
	float radius;
	float thickness;
	float ref_index;
	float aperture;
	float z_pos;
};

typedef struct lens_element lens_el;

typedef map<int, lens_el*> full_lens_data;

#endif
