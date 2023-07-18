// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ***************************************************************************
//
// file PMDiteration.cpp
//
// Purpose:
//			 PMDIteration class methods
//
// Programmer: Mathieu Lobet
// Creation:   Fri Oct 14 2016
//
// ***************************************************************************

#include "PMDIteration.h"

#include <DebugStream.h>

vector <PMDIteration> iterations;

// ***************************************************************************
// Method: PMDIteration::PMDIteration
//
// Purpose:
//			 Constructor
//
// Programmer: Mathieu Lobet
// Creation:   Fri Oct 14 2016
//
// Modifications:
//
// ***************************************************************************
PMDIteration::PMDIteration()
{
	this->dt		 = 0;
	this->time	     = 0;
	this->timeUnitSI = 1;
}

// ***************************************************************************
// Method: PMDIteration::~PMDIteration
//
// Purpose:
//			 Destructor
//
// Programmer: Mathieu Lobet
// Creation:   Fri Oct 14 2016
//
// Modifications:
//
// ***************************************************************************
PMDIteration::~PMDIteration()
{

}

// ***************************************************************************
// Method: PMDIteration::ScanFields
//
// Purpose:
//			 This method analyzes the group "fields" for the current iteration
//			 and gathers all the found datasets in the vector "fields".
//
// Arguments:
//			fileId : file id of the openPMD file used by hdf5
//
// Programmer: Mathieu Lobet
// Creation:   Fri Oct 14 2016
//
// Modifications:
//		 Mathieu Lobet, Tue Dec 13 2016
//		 I added the reading of the data Class and type size
// ***************************************************************************
void PMDIteration::ScanFields(hid_t fileId)
{
	hid_t				groupId;
	hid_t				subGroupId;
	hid_t				datasetId;
	hid_t				datasetType;
	string				path;
	hsize_t				nb_objects;
	hsize_t				nb_sub_objects;
	herr_t				err;
	int					i;
	int					i_sub_object;
	int					nb_attr;
	int					iattr;
	char				tmp_name[64];
	char				object_name[64];
	char				subObjectName[64];
	ssize_t				length;
	H5O_info_t			object_info;
	fieldGroupStruct	fieldGroup;
	PMDField			field;

	path = "/data/";
	// Add iteration to the global path
	path.append(name);
	path.append("/");
    // Add mesh path to the global path
	path.append(this->meshesPath);
	path.append("/");

	// Openning of the group "fields" of the current iteration
	groupId = H5Gopen2(fileId, path.c_str(), H5P_DEFAULT);


	// Get useful attributes fron the "fields" group
	field.ScanAttributes(groupId);

	/*
			Iteration over the objects in the group "fields"

			We look for the datasets and the groups.
			The idea is to make a vector containing all field datasets
			with their properties.
			I would prefer doing it recursively (and in a more beautiful manner)
			using the function H5Lvisit but the fact that useful attributes can
			be contained in groups make it more tricky and less readable.
			As a consequence, I prefer doing it by "hand".
			Be prepared for an explicit but ugly code...
			Mathieu
	*/

	  // First we get the number of objects
	  err = H5Gget_num_objs(groupId, &nb_objects);

	  // Iteration over the objects
	  for (i = 0; i < nb_objects; i++)
	  {
			// When an object is a group, we get useful information
			// from it and then scan the datasets inside.
			// If there are some sub-groups, they are not consider.

			// Get the object name
			length = H5Gget_objname_by_idx(groupId, (hsize_t)i,
				  object_name, (size_t) 64);

			// Get the type: group, dataset...
			// int object_type =  H5Gget_objtype_by_idx(groupId, (size_t)i );
			err = H5Oget_info_by_name(groupId, object_name , &object_info,
										H5P_DEFAULT);

			// Checking of the type
			switch(object_info.type)
			{
			// If group...
			case H5O_TYPE_GROUP:

				// Openning of the group
				subGroupId = H5Gopen2(groupId, object_name, H5P_DEFAULT);

				// Save the first part of the name
				strcpy(tmp_name,object_name);

				// Get useful attributes from the group
				field.ScanAttributes(subGroupId);

				// Create a new group structure
				strcpy(fieldGroup.name,object_name);
				fieldGroup.thetaComponents[0] = -1;
				fieldGroup.thetaComponents[1] = -1;
				fieldGroup.thetaComponents[2] = -1;
				fieldGroup.fieldIds.clear();

				// Copy the group geometry
				fieldGroup.geometry = field.geometry;

				  // Get the number of datasets
				  err = H5Gget_num_objs(subGroupId, &nb_sub_objects);

				// Then, we iterate over the datasets in this group
				for (i_sub_object = 0; i_sub_object < nb_sub_objects;
						 i_sub_object++)
				{
						// Get the dataset name
						length = H5Gget_objname_by_idx(subGroupId, (hsize_t)
													i_sub_object,
													subObjectName, (size_t) 64);

						// Save the name
						strcpy(field.name,tmp_name);
						strcat(field.name,"/");
						strcat(field.name,subObjectName);

						// Save the group path
						strcpy(field.groupPath,"/data/");
						strcat(field.groupPath,this->name.c_str());
						strcat(field.groupPath,"/");
						strcat(field.groupPath,this->meshesPath.c_str());
						strcat(field.groupPath,"/");
						strcat(field.groupPath,tmp_name);

						// Save the path
						strcpy(field.datasetPath,field.groupPath);
						strcat(field.datasetPath,"/");
						strcat(field.datasetPath,subObjectName);

						// Openning of the dataset
						datasetId = H5Dopen2(subGroupId,
											 subObjectName,
											 H5P_DEFAULT);

						// Get useful attributes from the group
						field.ScanAttributes(datasetId);

						// Set the grid dimension from the dataset
						field.SetGridDimensions(datasetId);

						// Get the class and size of data (float, integer...)
						datasetType = H5Dget_type(datasetId);
						field.dataSize = H5Tget_size(datasetType);
						field.dataClass = H5Tget_class(datasetType);

						// Insert the field in the list of files
						fields.push_back(field);

						// Store id of this field in field group
						fieldGroup.fieldIds.push_back(fields.size()-1);

						// Store id of known components for theta mode
						if (strcmp(subObjectName,"r")==0) // r component
						{
							fieldGroup.thetaComponents[0] = fields.size()-1;
						}
						else if (strcmp(subObjectName,"t")==0) //theta component
						{
							fieldGroup.thetaComponents[1] = fields.size()-1;
						}
						else if (strcmp(subObjectName,"z")==0) // z component
						{
							fieldGroup.thetaComponents[2] = fields.size()-1;
						}

						// Close the current dataSet
						H5Dclose(datasetId);

				  }

				  // Add the group to the list of field groups
				  fieldGroups.push_back(fieldGroup);

				  H5Gclose(subGroupId);

			break;
			// If dataset...
			case H5O_TYPE_DATASET:

				// Let's get the useful attributes of this dataset,
				// since it is localized in "fields"
				// it owns all its useful attributes.

				// Openning of the dataset
				datasetId = H5Dopen2(groupId, object_name, H5P_DEFAULT);

				// Save the name
				strcpy(field.name,object_name);

				// Save the group path
				strcpy(field.groupPath,"/data/");
				strcat(field.groupPath,this->name.c_str());
				strcat(field.groupPath,"/");
				strcat(field.groupPath,this->meshesPath.c_str());
				strcat(field.groupPath,"/");

				// Save the dataset path
				strcpy(field.datasetPath,field.groupPath);
				strcat(field.datasetPath,"/");
				strcat(field.datasetPath,object_name);

				// Scan useful attributes
				field.ScanAttributes(datasetId);

				// Set the grid dimension from the dataset
				field.SetGridDimensions(datasetId);

				// Get the class and size of data (float, integer...)
				datasetType = H5Dget_type(datasetId);
				field.dataSize = H5Tget_size(datasetType);
				field.dataClass = H5Tget_class(datasetType);

				// Insert the field in the vector of fields
				fields.push_back(field);

				// Close the current dataset
				H5Dclose(datasetId);

			break;

			default:
			debug5 << "visitLinks: node '" << this->name <<
			"' has an unknown type " << object_info.type << std::endl;
			break;
			}

	  }
	  H5Gclose(groupId);
}

// ***************************************************************************
// Method: PMDIteration::ScanParticles
//
// \brief This method analyzes the group "particles" for the current iteration.
//
//
// \author Programmer: Mathieu Lobet
// \date Creation:   Fri Oct 14 2016
//
// Modifications:
//
// ***************************************************************************
void PMDIteration::ScanParticles(hid_t fileId)
{

	int			i;
	string		path;
	char		objectName[64];
	hid_t		groupId;
	hid_t		particleGroupId;
	hsize_t		numObjects;
	ssize_t		length;
	herr_t		err = 1;
	H5O_info_t	objectInfo;

	// Creation of the path to the group "particles" having all the species
	path.append("/data/");
	path.append(this->name);
	path.append("/");
	path.append(this->particlesPath);
	path.append("/");

	if (this->particlesPath != "no_particles" && !this->particlesPath.empty())
	{
		err = H5Gget_objinfo (fileId, path.c_str(), 0, NULL);

		if (err!=0)
		{
			debug5 << " No group named particles of path: " << path << endl;
		}
		else
		{

			// Openning of the group "particles" of the current iteration
			groupId = H5Gopen2(fileId, path.c_str(), H5P_DEFAULT);

			// First we get the number of objects
			err = H5Gget_num_objs(groupId, &numObjects);

			// Iteration over the objects of the group "particles"
			for (i = 0; i < numObjects; i++)
			{

				// Create a temporary particle object
				PMDParticle			particle;

				// Get the particle group name
				length = H5Gget_objname_by_idx(groupId, (hsize_t)i,
				objectName, (size_t) 64);

				// Get info in order to get the type: group, dataset...
				err = H5Oget_info_by_name(groupId,
										objectName,
										&objectInfo,
										H5P_DEFAULT);

				// Check that the object is well a group
				if (H5O_TYPE_GROUP==objectInfo.type)
				{
					// Openning of the particle group
					particleGroupId = H5Gopen2(groupId,
											   objectName,
											   H5P_DEFAULT);

					// Save the name
					strcpy(particle.name,objectName);

					// Save the path
					strcpy(particle.path,path.c_str());
					strcat(particle.path,"/");
					strcat(particle.path,objectName);

					// Scan properties via datasets and attributes
					// in this particle group
					particle.ScanParticleGroup(particleGroupId);

					// Insert the particle object in the vector of particles
					this->particles.push_back(particle);

				}
			}
		}
	}
}

// ***************************************************************************
// Method: PMDIteration::HasFieldOfName
//
// Purpose:
//			 This method return true if the iteration object has a field
//			 component of specified name fieldName. This corresponds to
//			 the relative dataset name in the group `fields`.
//
// Inputs:
// fieldName: relative dataset name to be found
//
// Programmer: Mathieu Lobet
// Creation:   Thu Feb 09 2017
//
// Modifications:
//
// ***************************************************************************
bool PMDIteration::HasFieldOfName(char * fieldName)
{

  for (std::vector<PMDField>::iterator field = fields.begin() ;
		   field != fields.end(); ++field)
  {
	if (strcmp(field->name,fieldName))
	{
	  return true;
	}
  }

  return false;
}

template<typename T>
vector<T> PMDIteration::getAttributeArray(hid_t attrId, hid_t atype) {
	hid_t attrSpace = H5Aget_space(attrId);
	ssize_t arraySize = H5Sget_simple_extent_npoints(attrSpace);
	if (atype == H5T_STRING) {
		//if it is an array of strings muliply number of strings by string size
		atype = H5Aget_type(attrId);
		arraySize *= H5Tget_size(atype);
	}
	vector<T> attrArray(arraySize);
	H5Aread(attrId, atype, attrArray.data());
	return attrArray;
}

//TODO: maybe different name?
vector<string> PMDIteration::VectorCharToStr(vector<char> const& charVec,
											 size_t stringSize) {
	vector<string> stringVec;
	for (size_t i = 0; i < charVec.size(); i += stringSize) {
		stringVec.push_back(
			string(begin(charVec) + i, begin(charVec) + i + stringSize));
	}
	return stringVec;
}

//TODO: do something other than the "cout and return false" when encountering errors
//TODO: maybe move a lot of this into their own functions
bool PMDIteration::ReadAmrData(hid_t iterationId) {
	hid_t atype;
	hid_t attrId;

	//get number of Patches
	attrId = H5Aopen_name(iterationId, "numPatches");
	if (attrId >= 0) {
		long val = 0;
		H5Aread(attrId, H5T_NATIVE_LONG, &val); 
		patchChunks.resize(static_cast<size_t>(val));
	} else {cout << "not open numPatches"; return false;}

	//get the suffixes for all the patches
	vector<string> patchSuffixes;
	attrId = H5Aopen_name(iterationId, "patchSuffixes");
	if (attrId >= 0) {
		vector<char> charSuffixes = getAttributeArray<char>(attrId, H5T_STRING);
		size_t stringSize = H5Tget_size(H5Aget_type(attrId));
		patchSuffixes = VectorCharToStr(charSuffixes, stringSize);
	} else {cout << "not open patchSuffixes"; return false;}

	auto curPatchSuffix = begin(patchSuffixes);
	size_t domainNum = 0;
	for (vector<vector<chunk_t>>& patch : patchChunks) {
		vector<string> levelSuffixes;

		//get the suffixes for all the levels
		string levelSuffixName = "levelSuffixes" + *curPatchSuffix;
		attrId = H5Aopen_name(iterationId, levelSuffixName.c_str());
		if (attrId >= 0) {
			vector<char> charSuffixes = 
				getAttributeArray<char>(attrId, H5T_STRING);
			size_t stringSize = H5Tget_size(H5Aget_type(attrId));
			levelSuffixes = VectorCharToStr(charSuffixes, stringSize);
		} else {cout << "not open levelSuffixes"; return false;}

		//get all the chunk information from all the levels
		//assuming size of LevelSuffixes_patchxx == numLevels_patchxx
		for (string const& levelSuffix : levelSuffixes) {
			vector<chunk_t> chunks;
			string levelName = "chunkInfo" + levelSuffix;
			attrId = H5Aopen_name(iterationId, levelName.c_str());
			if (attrId >= 0) {
				vector<long> attrArray = getAttributeArray<long>(attrId, H5T_NATIVE_LONG);

				for (auto iter = begin(attrArray); iter != end(attrArray); iter += 6) {
					chunk_t newChunk;
					newChunk.lower[0] = iter[0];
					newChunk.lower[1] = iter[1];
					newChunk.lower[2] = iter[2];
					newChunk.upper[0] = iter[3];	
					newChunk.upper[1] = iter[4];
					newChunk.upper[2] = iter[5];
					newChunk.domainNumber = domainNum;
					domainNum++;
					chunks.push_back(newChunk);
					// cout << newChunk << endl;
				}
			} else {cout << "not open chunkInfo"; return false;}
			patch.push_back(chunks);
		}
		curPatchSuffix++;
	}

	ReadChildList();

	return true;
}

void PMDIteration::ReadChildList() {
	for (vector<vector<chunk_t>>& patch : patchChunks) {
		for (size_t level = 0; level < patch.size() - 1; ++level) {
			for (chunk_t& chunk : patch[level]) {
				chunk.childList =
					FindChildListForChunk(chunk, patch[level + 1]);
				cout << chunk << endl;
			}
		}
	}
}


//TODO: currently assumes that chunks in "level" are more fine than "chunk"
vector<size_t> PMDIteration::FindChildListForChunk(const chunk_t& chunk,
		const vector<chunk_t>& level) const {
	vector<size_t> childList;
	for (const chunk_t& levelChunk : level) {
		bool isInChunk = true;
		for (size_t side = 0; side < 3; ++side) {
			if (levelChunk.upper[side] < chunk.lower[side] || 
					levelChunk.lower[side] > chunk.upper[side]) {
				isInChunk = false;
				break;
			}
		}
		if (isInChunk) {
			childList.push_back(levelChunk.domainNumber);
		}
	}
	return childList;
}

vector<chunk_t> PMDIteration::getChunk(size_t patchNum, size_t levelNum) const {
	cout << patchChunks.size() << endl;
	cout << patchChunks[patchNum].size() << endl;
	cout << patchChunks[patchNum][levelNum].size() << endl;
	return patchChunks[patchNum][levelNum];
}

// ***************************************************************************
// Method: PMDIteration::PrintInfo
//
// Purpose:
//			 This method prints the content of an object PMDIteration:
//			 - field datasets
//			 - particle datasets
//
// Programmer: Mathieu Lobet
// Creation:   Fri Oct 14 2016
//
// Modifications:
//
// ***************************************************************************
void PMDIteration::PrintInfo()
{

	int i;

	cout << " Information about iteration " << name << endl;
	cout << " - dt: "<< dt << endl;
	cout << " - time: "<< time << endl;
	cout << " - timeUnitSI: "<< timeUnitSI << endl;
	cout << endl;
	cout << " Number of field datasets: " << fields.size() << endl;
	for (std::vector<PMDField>::iterator field = fields.begin() ;
			 field != fields.end(); ++field)
	{
		  cout << " Field name: " << field->name << endl;
		  cout << " - dataset path: " << field->datasetPath << endl;
		  cout << " - group path: " << field->groupPath << endl;
		  cout << " - dx: " << field->gridSpacing[0]
			   << " dy: " << field->gridSpacing[1]
			   << " dz: " << field->gridSpacing[2] << endl;
		  cout << " - xmin: " << field->gridGlobalOffset[0]
			   << " ymin: " << field->gridGlobalOffset[1]
			   << " zmin: " << field->gridGlobalOffset[2] << endl;
		  cout << " - number of dimensions: "
			   << field->ndims << endl;
		  cout << " - nx: " << field->nbNodes[0]
			   << " ny: " << field->nbNodes[1]
			   << " nz: " << field->nbNodes[2] << endl;
		  cout << " - xshift: " << field->gridPosition[0]
			   << " yshift: " << field->gridPosition[1]
			   << " zshift: " << field->gridPosition[2] << endl;
		  cout << " - Unit SI: " << field->unitSI<< endl;
		  cout << " - Grid Unit SI: " << field->gridUnitSI<< endl;
		  cout << " - Geometry: " << field->geometry.c_str() << endl;
		  cout << " - xlabel: " << field->axisLabels[0].c_str()
			   << " ylabel: " << field->axisLabels[1].c_str()
			   << " zlabel: " << field->axisLabels[2].c_str() << endl;
		  cout << " - Units: " << field->unitsLabel.c_str() << endl;
		  cout << " - Data order: " << field->dataOrder.c_str() << endl;
		  cout << endl;
	}
	cout << endl;
	  cout << " Number of particle groups: " << particles.size() << endl;
	  for (std::vector<PMDParticle>::iterator particle = particles.begin() ;
			 particle != particles.end(); ++particle)
	{
		  cout << endl;
		  cout << " Particle name: " << particle->name << endl;
		  cout << " - charge: " << particle->charge << endl;
		  cout << " - mass: " << particle->mass << endl;
		  cout << " - number of particles: " << particle->numParticles << endl;
		  cout << " - Position datasets: " << particle->positionsId[0]
			   << " " << particle->positionsId[1]
			   << " "<< particle->positionsId[2]
			   << endl;
		  cout << " Number of scalar datasets: "
			   << particle->GetNumScalarDatasets() << endl;
			for (i=0;i<particle->GetNumScalarDatasets();i++)
		  {
				cout << " - Name: " << particle->scalarDataSets[i].name
                    << ", Unit Label: " << particle->scalarDataSets[i].unitLabel.c_str()
					<< ", path: " << particle->scalarDataSets[i].path
					<< ", unitSI: " << particle->scalarDataSets[i].unitSI
					<< endl;
		  }
		  cout << " Number of vector datasets: "
			   << particle->GetNumVectorDatasets() << endl;
			for (i=0;i<particle->GetNumVectorDatasets();i++)
		  {
				cout << " - Name: " << particle->vectorDataSets[i].name
                     << ", Unit Label: " << particle->vectorDataSets[i].unitLabel.c_str()
					 << ", path: " << particle->vectorDataSets[i].path
					 << ", scalar datasets: "
					 << particle->vectorDataSets[i].dataSetId[0]
					 << " " << particle->vectorDataSets[i].dataSetId[1]
					 << " " << particle->vectorDataSets[i].dataSetId[2]
					 << endl;
		  }
	}
}
