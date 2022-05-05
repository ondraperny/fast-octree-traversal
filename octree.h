// ===================================================================
//
// octree.h
//      OCTREE - Octree implementation with an efficient parametric algorithm for octree traversal
//
// Class: OCTREE
//

#ifndef __OCTREE__
#define __OCTREE__

// standard C++ headers
#include <fstream>
#include <chrono>

// forward declarations
class CObjectList;
class CObject3D;
class CRay;

// ----------------------------------------------------------
/// Structure for building and querying octree with fast parametric algorithm
class OCTREE :
	public CASDS_BB
{
public:
	OCTREE();
	virtual ~OCTREE();

protected:
	/// Octree node - terminal node holds pointer to vector, non-terminal node hold pointer to array of 8 pointers to children nodes
	struct OctreeNode {
		bool terminal = true;
		union {
			OctreeNode** childrenNodes;
			std::vector<CObject3D*>* objects;
		};
	};

	/// Struct for informations needed for parametric algorithm using queue
	struct queueNode {
		queueNode(float tx0, float ty0, float tz0, float tx1, float ty1, float tz1, OctreeNode* node) {
			this->tx0 = tx0;
			this->ty0 = ty0;
			this->tz0 = tz0;
			this->tx1 = tx1;
			this->ty1 = ty1;
			this->tz1 = tz1;
			this->node = node;
		}
		float tx0;
		float ty0;
		float tz0;
		float tx1;
		float ty1;
		float tz1;
		OctreeNode* node;
	};
	// PARAMETERS
	const int maxDepth = 9;

	// Functions BASE
	OctreeNode* root = nullptr;
	
	/// build the Octree structure
	virtual void BuildUp(const CObjectList& objlist);
	virtual void GetBBox(SBBox& box) { box = bbox; }
	/// return ID
	virtual void ProvideID(ostream& app);
	/// makes necessary call when destroying octree
	virtual void Remove();
	/// return Octree ID
	virtual EASDS_ID GetID() const { return ID_OctreeTraversal; }
	/// recursively build the tree
	OctreeNode* BuildOctreeRecursively(const std::vector<CObject3D*> objlist,
		const int currentOctreeDepth,
		const SBBox& currentBounds);
	/// Recursively delete whole Octree (traversing nodes)
	void deleteTreeRec(OctreeNode* current);
	/// find nearest intersection
	virtual const CObject3D* FindNearestI(CRay& ray, CHitPointInfo& info);

	// parametric ALGORITHM functions and variables
	unsigned char a = 0;
	CRay* currentRay = nullptr;
	CHitPointInfo* currentInfo = nullptr;
	CObject3D* finalObj = nullptr;
	bool stopRay = false;

	/// Calculate which sub-node is first intersected - therefore 
	/// the sub-node which ray use to enter the parent node
	inline int first_node(float tx0, float ty0, float tz0, float txm, float tym, float tzm);
	/// Based on input parameters calculate next intersected node
	inline int new_node(float tx1, int x, float ty1, int y, float tz1, int z);
	/// Correct direction of ray (if any component is negative) and calculate parameters for root node
	void ray_parameter(OctreeNode* node, CRay ray);
	/// Main function processing the nodes (recursive)
	void proc_subtree(float tx0, float ty0, float tz0, float tx1, float ty1, float tz1, OctreeNode* node);
	/// Main function processing the nodes (non-recursive)
	void proc_subtreeNonRec(float tx0, float ty0, float tz0, float tx1, float ty1, float tz1, OctreeNode* node);

	// Functions STATISTICS
	std::chrono::high_resolution_clock::time_point clock;
	unsigned int numberIntersectTest = 0;
	unsigned int numberOfQueries = 0;
	unsigned int numberOfTraverStep = 0;
	std::chrono::milliseconds buildTime;

	// Functions HELPER/DEBUG
	ofstream myfile;
	unsigned int numberOfLeaves = 0;
	bool endSearch = false;
	unsigned int termNodesCnt = 0;

	// Functions for GLV visualization
	/// Print BBox as rectangle (with specified color and lineWidth) to the file for GLV visualization
	void printVerticesFromBBox(const SBBox& boxPrint, string RGBcolor = "0.6 0.6 0.6", string lineWidth = "1");
	/// Write triangle (with specified color) to the file for GLV visualization
	void printTriangle(CObject3D& triangle, string color = "0 1 0");
	/// Write ray to the file for GLV visualization
	void printRay(CRay& ray);
	/// Write all triangles to the file for GLV visualization
	void printAllMeshes(const CObjectList& objlist);
	/// Write all triangles (with specified color) to the file for GLV visualization
	void printAllMeshes(std::vector<CObject3D*> objlist, string color = "1 1 1");

	// Functions for debug console output
	/// Traverse octree for debugging and printing purposes
	void traverseOctree(OctreeNode* current);
	/// Print to console information about object position
	void printTriagle(CObject3D* object);
	/// Print to console information of all objects in CObjectList for debugging purposes
	void printObjects(const CObjectList& objlist);
};

#endif // __OCTREE__

