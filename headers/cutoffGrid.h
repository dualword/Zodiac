/*   Zodiac - molecular modelling environment. www.zeden.org
    Copyright (C) 2008  Nicola Zonta (nicola.zonta(at)zeden.org)

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

#ifndef cutoffGrid_H
#define cutoffGrid_H

#include <vector>

using namespace std;

template <class T>
struct objectList {
	vector<T> objects;
};

template <class T>
struct staticCutoffNode {
	vector<T> cellObjects;
	objectList<T>* objects;
};

template <class T>
class cutoffGrid
{
public:
  	cutoffGrid(vector<T>& objects, float cutRadius);
    ~cutoffGrid();
    objectList<T>* getNeighborObjects(vect pos);
    vector<staticCutoffNode<T>*> grid;
    int gridSize[3];
	float gridResolution;
	float gridMin[3];
	int gridSizeXY;
	
private:
	float cutoffRadius;
	float gridMax[3];
	unsigned int maxSID;
	
};


template <class T>
cutoffGrid<T>::cutoffGrid(vector<T>& objects, float cutRadius)
{

	gridResolution = cutRadius;
	
	cutoffRadius = cutRadius;
	
	if (!objects.size()) {
		gridMin[0] = 0.0;	gridMin[1] = 0.0;	gridMin[2] = 0.0;
		gridSize[0] = 0;	gridSize[1] = 0;	gridSize[2] = 0;
	}
	else {
	
		for (unsigned int i = 0; i < objects.size(); i++) {
			if (!i) {
				gridMin[0] = objects[i]->GetVector ().x();
				gridMin[1] = objects[i]->GetVector ().y();
				gridMin[2] = objects[i]->GetVector ().z();

				gridMax[0] = objects[i]->GetVector ().x();
				gridMax[1] = objects[i]->GetVector ().y();
				gridMax[2] = objects[i]->GetVector ().z();
			}
			else {
				if (objects[i]->GetVector ().x() < gridMin[0]) {
					gridMin[0] = objects[i]->GetVector ().x();
				}
				if (objects[i]->GetVector ().y() < gridMin[1]) {
					gridMin[1] = objects[i]->GetVector ().y();
				}
				if (objects[i]->GetVector ().z() < gridMin[2]) {
					gridMin[2] = objects[i]->GetVector ().z();
				}
				if (objects[i]->GetVector ().x() > gridMax[0]) {
					gridMax[0] = objects[i]->GetVector ().x();
				}
				if (objects[i]->GetVector ().y() > gridMax[1]) {
					gridMax[1] = objects[i]->GetVector ().y();
				}
				if (objects[i]->GetVector ().z() > gridMax[2]) {
					gridMax[2] = objects[i]->GetVector ().z();
				}
				
			}
		}

		gridMin[0] -= 2.0 * cutRadius;
		gridMin[1] -= 2.0 * cutRadius;
		gridMin[2] -= 2.0 * cutRadius;	  
		
		gridMax[0] += 2.0 * cutRadius;
		gridMax[1] += 2.0 * cutRadius;
		gridMax[2] += 2.0 * cutRadius;	  


		gridSize[0] = (int)((gridMax[0] - gridMin[0]) / gridResolution + 0.5);
		gridSize[1] = (int)((gridMax[1] - gridMin[1]) / gridResolution + 0.5);
		gridSize[2] = (int)((gridMax[2] - gridMin[2]) / gridResolution + 0.5);
		
		

		gridSizeXY = gridSize[0] * gridSize[1];

		grid.resize(gridSize[0] * gridSize[1] * gridSize[2], NULL);

		vector<int*> nodes;

		for (unsigned int i = 0; i < objects.size(); i++) {

			int gridPos[3];

			gridPos[0] = (int)((objects[i]->GetVector ().x() - gridMin[0]) / gridResolution + 0.5);
			gridPos[1] = (int)((objects[i]->GetVector ().y() - gridMin[1]) / gridResolution + 0.5);
			gridPos[2] = (int)((objects[i]->GetVector ().z() - gridMin[2]) / gridResolution + 0.5);

			if (gridPos[0] >= 0 && gridPos[0] < gridSize[0] &&
				gridPos[1] >= 0 && gridPos[1] < gridSize[1] &&
				gridPos[2] >= 0 && gridPos[2] < gridSize[2]) {

				int index = gridPos[0] + gridPos[1] * gridSize[0] + gridPos[2] * gridSizeXY;
				if (grid[index] == NULL) {
					grid[index] = new staticCutoffNode<T>;
					grid[index]->objects = new objectList<T>;
					grid[index]->cellObjects.push_back(objects[i]);
					int *nodeIndex = new int[3];
					nodeIndex[0] = gridPos[0];
					nodeIndex[1] = gridPos[1];
					nodeIndex[2] = gridPos[2];

					nodes.push_back(nodeIndex);
				}
				else {	
					grid[index]->cellObjects.push_back(objects[i]);
				}
			}
		}

		int radius = 1;

		int x, y, z;

		vector<int*> sphere;

		for (x = -radius; x <= radius; x++)
			for (y = -radius; y <= radius; y++)
				for (z = -radius; z <= radius; z++) {
					int *p = new int[3];			
					p[0] = x;
					p[1] = y;
					p[2] = z;
					sphere.push_back(p);
				}
		for (unsigned int i = 0; i < nodes.size(); i++) {
			int *ni = nodes[i];
			int index = ni[0] + ni[1] * gridSize[0] + ni[2] * gridSizeXY;
			for (unsigned int j = 0; j < sphere.size(); j++) {
				int *p = sphere[j];
				int newPoint[3];
				newPoint[0] = ni[0] + p[0];
				newPoint[1] = ni[1] + p[1];
				newPoint[2] = ni[2] + p[2];

				if (newPoint[0] >= 0 && newPoint[0] < gridSize[0] &&
					newPoint[1] >= 0 && newPoint[1] < gridSize[1] &&
					newPoint[2] >= 0 && newPoint[2] < gridSize[2]) {
					int newIndex = newPoint[0] + newPoint[1] * gridSize[0] + newPoint[2] * gridSizeXY;
					if (grid[newIndex] == NULL) {
						grid[newIndex] = new staticCutoffNode<T>;
						grid[newIndex]->objects = new objectList<T>;
					}
					for (unsigned int k = 0; k < grid[index]->cellObjects.size(); k++) {
						grid[newIndex]->objects->objects.push_back(grid[index]->cellObjects[k]);
					}
				}
			}	
			delete[] nodes[i];
		}
		
		for (unsigned int i = 0; i < sphere.size(); i++) {
			if (sphere[i]) {
				delete[] sphere[i];
				sphere[i] = NULL;
			}
		}
		sphere.clear();
	}

}


template <class T>
cutoffGrid<T>::~cutoffGrid()
{
	for (unsigned int i = 0; i < grid.size(); i++) {
		if (grid[i]) {
			if (grid[i]->objects) {
				grid[i]->objects->objects.clear();
				delete grid[i]->objects;
			}
			grid[i]->objects = NULL;
			delete grid[i];
			grid[i] = NULL;
		}
	}
	grid.clear();
}

template <class T>
objectList<T>* cutoffGrid<T>::getNeighborObjects(vect pos)
{
	int gridPos[3];
	
	gridPos[0] = (int)((pos.x() - gridMin[0]) / gridResolution + 0.5);
	gridPos[1] = (int)((pos.y() - gridMin[1]) / gridResolution + 0.5);
	gridPos[2] = (int)((pos.z() - gridMin[2]) / gridResolution + 0.5);
	
	if (gridPos[0] >= 0 && gridPos[0] < gridSize[0] &&
		gridPos[1] >= 0 && gridPos[1] < gridSize[1] &&
		gridPos[2] >= 0 && gridPos[2] < gridSize[2]) {
		
		int index = gridPos[0] + gridPos[1] * gridSize[0] + gridPos[2] * gridSizeXY;
		if (grid[index]) {
			return grid[index]->objects;
		}
		else {
			return NULL;
		}
	}
	return NULL;
}

#endif

