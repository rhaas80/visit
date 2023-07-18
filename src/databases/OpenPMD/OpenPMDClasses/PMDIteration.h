// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef PMDITERATION_H
#define PMDITERATION_H

#include "PMDField.h"
#include "PMDParticle.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#include <hdf5.h>
//#include <visit-hdf5.h>

//#include "H5Cpp.h"

using namespace std;

struct chunk_t {
	long lower[3]; //lower end of chunk
	long upper[3]; //upper end of chunk
};

// ***************************************************************************
// Class: PMDIteration
//
// Purpose:
//	  This class enables to manage the different iterations of
//	  an openPMD file.
//	  When an openPMD file is read, a PMDIteration is created
//	  for each iteration.
//
// Programmer: Mathieu Lobet
// Creation:   Fri Oct 14 2016
//
// Modifications:
//
// ***************************************************************************
class PMDIteration
{
	public:
	PMDIteration();
	~PMDIteration();

	void	ScanFields(hid_t fileId);
	void	ScanParticles(hid_t fileId);
	void	PrintInfo();
	bool	HasFieldOfName(char * fieldName);
	bool	ReadAmrData(hid_t iterationId);
	const vector<chunk_t>& getChunk(size_t patchNum, size_t levelNum) const;

	// Iteration attributes
	/// Iteration name
	string	name;
	/// Mesh path in the iteration group
	string	meshesPath;
	/// Particles path in the iteration group
	string	particlesPath;
	/// Iteration time step
	float  	dt;
	/// Iteration corresponding time
	float 	time;
	/// factor to convert the time in SI units
	float 	timeUnitSI;

	/// Vector of field objects from the datasets
	vector <PMDField> fields;

	/// Vector of field group structures
	vector <fieldGroupStruct> fieldGroups;

	/// Vector of particle objects
	vector <PMDParticle> particles;

	protected:

	private:
	vector<vector<vector<chunk_t>>> patchChunks;

	template<typename T>
	vector<T> getAttributeArray(hid_t attrId, hid_t atype);
	vector<string> VectorCharToStr(vector<char> const& charVec,
								   size_t stringSize);
};

inline std::ostream& operator<<(std::ostream& output, chunk_t const& chunk) {
	output << "lower:" << endl;
	output << "\t" << chunk.lower[0] << ", " 
		   		   << chunk.lower[1] << ", "
				   << chunk.lower[2] << endl;
	output << "upper:" << endl;
	output << "\t" << chunk.upper[0] << ", "
				   << chunk.upper[1] << ", "
				   << chunk.upper[2];
	return output;
}

#endif
