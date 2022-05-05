// standard C++ headers
#include <iomanip>
#include <queue>

#include "octree.h"

#define GENERATE_MODEL false
#define DEBUG_PRINT false

// ----------------------------------------------------------
// class OCTREE

// default constructor
OCTREE::OCTREE()
{}

// default destructor
OCTREE::~OCTREE()
{}

void OCTREE::Remove()
{
    std::cout.precision(5);
    // print statistics before removing tree
    std::cout << "Complete time (build and traversal):    " << world->totalTime + buildTime.count() << endl;
    std::cout << "Traversal time:                         " << world->totalTime << endl;
    std::cout << "T_B (structure build time):             " << buildTime.count() << endl;
    std::cout << "T_R (average time per query):           " << world->queriesTime << endl;
    std::cout << "N_IT (average incidence per query):     " << (double)numberIntersectTest / (double)numberOfQueries << endl;
    std::cout << "N_TR (average travers steps per query): " << (double)numberOfTraverStep / (double)numberOfQueries << endl;
    std::cout << "N_Q (number of queries):                " << numberOfQueries << std::endl;
    
    if (GENERATE_MODEL)
        myfile.close();

    deleteTreeRec(root);
}

void OCTREE::deleteTreeRec(OctreeNode* currentNode) {
    if (!currentNode->terminal) {
        for (size_t i = 0; i < 8; i++) {
            if (currentNode->childrenNodes[i]) {
                deleteTreeRec(currentNode->childrenNodes[i]);
            }
        }
    }
    if (currentNode->terminal)
        delete currentNode->objects;
    else 
        delete[] currentNode->childrenNodes;
    delete currentNode;
}

void OCTREE::BuildUp(const CObjectList& objlist)
{
    // Initialize the box of the whole scene
    bbox.Initialize();
    std::cout << "Box init: " << bbox.Min() << bbox.Max() << std::endl;

    // Initial box of the whole scene
    InitializeBox(this->bbox, (CObjectList&)objlist);
    std::cout << "Box init: " << bbox.Min() << bbox.Max() << std::endl;

    CVector3D boxSize = bbox.Max() - bbox.Min();
    CVector3D centroid = bbox.ComputeCentroid();
    float maxSize = max(max(abs(boxSize.x), abs(boxSize.y)), abs(boxSize.z));

    SBBox firstBox(CVector3D(0, 0, 0), CVector3D(abs(bbox.MaxX()), abs(bbox.MaxY()), abs(bbox.MaxZ())));
    SBBox seconBox(CVector3D(0, 0, 0), CVector3D(abs(bbox.MinX()), abs(bbox.MinY()), abs(bbox.MinZ())));
    float edgeLen = max(firstBox.Diagonal().MaxComponent(), seconBox.Diagonal().MaxComponent());

    CVector3D newDiagonal(edgeLen, edgeLen, edgeLen);
    CVector3D halfDiagonal = CVector3D(maxSize / 2, maxSize / 2, maxSize / 2);

    bbox.SetMin(CVector3D(0, 0, 0) - newDiagonal);
    bbox.SetMax(CVector3D(0, 0, 0) + newDiagonal);
    
    // start time cound for building octree
    clock = std::chrono::high_resolution_clock::now();
    // build octree
    root = BuildOctreeRecursively(objlist, 0, bbox);

    // print time needed for build
    //std::cout << "T_B (structure build time): " << (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clock)).count() <<
    //    " milliseconds" << std::endl;
    buildTime = (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clock));

    // traversing and showing the tree
    if (GENERATE_MODEL) {
        traverseOctree(root);
        myfile.open("example.txt");
    }

    std::cout << "Number of objects:" << objlist.size() << std::endl;
    std::cout << "Number of leaves:" << numberOfLeaves << std::endl;

    STATUS << "OCTREE construction finished" << endl;
    clock = std::chrono::high_resolution_clock::now();

    return;
}

OCTREE::OctreeNode* OCTREE::BuildOctreeRecursively(const std::vector<CObject3D*> objlist,
    const int currentOctreeDepth,
    const SBBox& currentBounds) {
    OctreeNode* currentNode = new OctreeNode();

    std::vector<CObject3D*> childObjects[8];
    unsigned int childrenObjectsNumber[8] = { 0 }; // initialize to 0s

    /// print output
    if (DEBUG_PRINT) {
        CVector3D tmp = currentBounds.Max() - currentBounds.Min();
        std::cout << "Min/Max: " << currentBounds.Min() << " " << currentBounds.Max() << std::endl;
        std::cout << "Distance: " << sqrtf(powf(tmp.getX(), 2) + powf(tmp.getY(), 2) + powf(tmp.getZ(), 2)) << std::endl;
    }

    CVector3D center;
    currentBounds.ComputeCentroid(center);
    CVector3D sideLengths = currentBounds.Max() - currentBounds.Min();

    // calculate positions of children nodes for intersection test
    SBBox childBBoxes[8] = {
    SBBox(currentBounds.Min(), center),
    SBBox((currentBounds.Min() + CVector3D(0, 0, sideLengths.z / 2.0)), center + CVector3D(0, 0, sideLengths.z / 2.0)), 
    SBBox((currentBounds.Min() + CVector3D(0, sideLengths.y / 2.0, 0)), center + CVector3D(0, sideLengths.y / 2.0, 0)),
    SBBox((center - CVector3D(sideLengths.x / 2.0, 0, 0)), currentBounds.Max() - CVector3D(sideLengths.x / 2.0, 0, 0)),
    SBBox((currentBounds.Min() + CVector3D(sideLengths.x / 2.0, 0, 0)), center + CVector3D(sideLengths.x / 2.0, 0, 0)),
    SBBox((center - CVector3D(0, sideLengths.y / 2.0, 0)), currentBounds.Max() - CVector3D(0, sideLengths.y / 2.0, 0)),
    SBBox((center - CVector3D(0, 0, sideLengths.z / 2.0)), currentBounds.Max() - CVector3D(0, 0, sideLengths.z / 2.0)),
    SBBox(center, currentBounds.Max()), 
    };

    // iterate objects in scene and assign them to corresponding nodes
    for (auto& object : objlist)
    {
        for (int i = 0; i < 8; i++)
        {
            if (OverlapS(object->GetBox(), childBBoxes[i])) {
                childObjects[i].push_back(object);
                childrenObjectsNumber[i]++;
            }
        }
    }

    if (currentOctreeDepth < maxDepth/* && objlist.size() >= 10*/) {
        for (auto& childObjNumber : childrenObjectsNumber) {
            if (childObjNumber > 0) {
                currentNode->terminal = false;
                break;
            }
        }

        // allocate array for children nodes
        currentNode->childrenNodes = new OctreeNode * [8]{ nullptr };

        for (int i = 0; i < 8; i++)
        {
            if (childrenObjectsNumber[i] > 0) {
                //currentNode->childrenCount += 1;
                //currentNode->terminal = false;
                currentNode->childrenNodes[i] = BuildOctreeRecursively(childObjects[i], currentOctreeDepth + 1, childBBoxes[i]);
            }
            else
            {
                currentNode->childrenNodes[i] = nullptr;
                if (DEBUG_PRINT)
                    std::cout << "No children " << i << ", depth: " << currentOctreeDepth << std::endl;
            }
        }
    }
    else {
        // allocate vector for objects
        currentNode->objects = new std::vector<CObject3D*>();
        *(currentNode->objects) = objlist;
        currentNode->objects->shrink_to_fit();

        if (DEBUG_PRINT)
            std::cout << "Number of objects in current node: " << currentNode->objects->size() << std::endl;
        numberOfLeaves += 1;
    }
    //std::cout << sizeof(*currentNode) << std::endl;
    return currentNode;
}

inline int OCTREE::first_node(float tx0, float ty0, float tz0, float txm, float tym, float tzm) {
    // 8 bites to remember from which side ray intersected parent node - choosing first sub-node
    unsigned char bitPosition = 0;

    // select the entry plane and set bits
    // tx0 > ty0 - YZ plane, possible affected bits 1, 0
    if (tx0 > ty0) {
        if (tx0 > tz0) {
            if (tym < tx0)
                bitPosition |= 2; 
            if (tzm < tx0)
                bitPosition |= 1;
            return static_cast<int>(bitPosition);
        }
    }
    // (tx0 < ty0) - XZ plane, possible affected bits 2, 0
    else {
        if (ty0 > tz0) { 
            if (txm < ty0)
                bitPosition |= 4;
            if (tzm < ty0)
                bitPosition |= 1;
            return static_cast<int>(bitPosition);
        }
    }
    // if tx0 < ty0 and ty0 < tz0 - XY plane, possible affected bits 2, 1
    if (txm < tz0)
        bitPosition |= 4;
    if (tym < tz0)
        bitPosition |= 2;
    return static_cast<int>(bitPosition);
}

inline int OCTREE::new_node(float tx1, int x, float ty1, int y, float tz1, int z) {
    // choosing the entry plane
    // if txm is smallest parameter, then entry plane is YZ, therefore return x, same for rest
    if (tx1 < ty1) {
        if (tx1 < tz1)
            return x;
    }
    else {
        if (ty1 < tz1)
            return y;
    }
    return z;
}


void OCTREE::proc_subtree(float tx0, float ty0, float tz0, float tx1, float ty1, float tz1, OctreeNode* node) {
    float txm, tym, tzm;
    int currNode;

    // if algorithm try to go the node that does exist in octree - space of that potential node is empty
    // or if any t < 0 -> ray miss the node
    if (node == nullptr || tx1 < 0 || ty1 < 0 || tz1 < 0 || termNodesCnt > 3) return;
    numberOfTraverStep++;

    // if reach terminal node, check if there is intersection between ray and object 
    if (node->terminal) {
        for (auto& object : *(node->objects)) {
            numberIntersectTest++;
            if (object->NearestInt(*currentRay, *currentInfo)) {
                currentInfo->maxt = currentInfo->t;
                finalObj = object;
            }
        }
        if (finalObj) {
            endSearch = true;
            termNodesCnt++;
        }
        return;
    }

    txm = 0.5 * (tx0 + tx1);
    tym = 0.5 * (ty0 + ty1);
    tzm = 0.5 * (tz0 + tz1);

    currNode = first_node(tx0, ty0, tz0, txm, tym, tzm);
    do {
        switch (currNode)
        {
        case 0: {
            proc_subtree(tx0, ty0, tz0, txm, tym, tzm, node->childrenNodes[a]);
            currNode = new_node(txm, 4, tym, 2, tzm, 1);
            break; }
        case 1: {
            proc_subtree(tx0, ty0, tzm, txm, tym, tz1, node->childrenNodes[1^a]);
            currNode = new_node(txm, 5, tym, 3, tz1, 8);
            break; }
        case 2: {
            proc_subtree(tx0, tym, tz0, txm, ty1, tzm, node->childrenNodes[2^a]);
            currNode = new_node(txm, 6, ty1, 8, tzm, 3);
            break; }
        case 3: {
            proc_subtree(tx0, tym, tzm, txm, ty1, tz1, node->childrenNodes[3^a]);
            currNode = new_node(txm, 7, ty1, 8, tz1, 8);
            break; }
        case 4: {
            proc_subtree(txm, ty0, tz0, tx1, tym, tzm, node->childrenNodes[4^a]);
            currNode = new_node(tx1, 8, tym, 6, tzm, 5);
            break; }
        case 5: {
            proc_subtree(txm, ty0, tzm, tx1, tym, tz1, node->childrenNodes[5^a]);
            currNode = new_node(tx1, 8, tym, 7, tz1, 8);
            break; }
        case 6: {
            proc_subtree(txm, tym, tz0, tx1, ty1, tzm, node->childrenNodes[6^a]);
            currNode = new_node(tx1, 8, ty1, 8, tzm, 7);
            break; }
        case 7: {
            proc_subtree(txm, tym, tzm, tx1, ty1, tz1, node->childrenNodes[7^a]);
            currNode = 8;
            break; }
        }
    } while (currNode < 8);
}

void OCTREE::proc_subtreeNonRec(float tx0, float ty0, float tz0, float tx1, float ty1, float tz1, OctreeNode* node) {
    float txm, tym, tzm;
    int currNode;

    std::queue<queueNode> nodeQueue;
    queueNode root(tx0, ty0, tz0, tx1, ty1, tz1, node);
    nodeQueue.push(root);

    while (!nodeQueue.empty() && termNodesCnt <= 3) {
        tx0 = nodeQueue.front().tx0;
        ty0 = nodeQueue.front().ty0;
        tz0 = nodeQueue.front().tz0;
        tx1 = nodeQueue.front().tx1;
        ty1 = nodeQueue.front().ty1;
        tz1 = nodeQueue.front().tz1;
        node = nodeQueue.front().node;
        nodeQueue.pop();

        // if algorithm try to go the node that does exist in octree - space of that potential node is empty
        // or if any t < 0 -> ray miss the node
        if (node == nullptr || tx1 < 0 || ty1 < 0 || tz1 < 0)
            continue;

        numberOfTraverStep++;
        // if reach terminal node, check if there is intersection between ray and object 
        if (node->terminal) {
            for (auto& object : *(node->objects)) {
                numberIntersectTest++;
                if (object->NearestInt(*currentRay, *currentInfo)) {
                    currentInfo->maxt = currentInfo->t;
                    finalObj = object;
                }
            }
            if (finalObj) {
                endSearch = true;
                termNodesCnt++;
            }
            continue;
        }

        txm = 0.5 * (tx0 + tx1);
        tym = 0.5 * (ty0 + ty1);
        tzm = 0.5 * (tz0 + tz1);

        currNode = first_node(tx0, ty0, tz0, txm, tym, tzm);
        do {
            switch (currNode)
            {
            case 0: {
                nodeQueue.push(queueNode(tx0, ty0, tz0, txm, tym, tzm, node->childrenNodes[a]));
                currNode = new_node(txm, 4, tym, 2, tzm, 1);
                break; }
            case 1: {
                nodeQueue.push(queueNode(tx0, ty0, tzm, txm, tym, tz1, node->childrenNodes[1 ^ a]));
                currNode = new_node(txm, 5, tym, 3, tz1, 8);
                break; }
            case 2: {
                nodeQueue.push(queueNode(tx0, tym, tz0, txm, ty1, tzm, node->childrenNodes[2 ^ a]));
                currNode = new_node(txm, 6, ty1, 8, tzm, 3);
                break; }
            case 3: {
                nodeQueue.push(queueNode(tx0, tym, tzm, txm, ty1, tz1, node->childrenNodes[3 ^ a]));
                currNode = new_node(txm, 7, ty1, 8, tz1, 8);
                break; }
            case 4: {
                nodeQueue.push(queueNode(txm, ty0, tz0, tx1, tym, tzm, node->childrenNodes[4 ^ a]));
                currNode = new_node(tx1, 8, tym, 6, tzm, 5);
                break; }
            case 5: {
                nodeQueue.push(queueNode(txm, ty0, tzm, tx1, tym, tz1, node->childrenNodes[5 ^ a]));
                currNode = new_node(tx1, 8, tym, 7, tz1, 8);
                break; }
            case 6: {
                nodeQueue.push(queueNode(txm, tym, tz0, tx1, ty1, tzm, node->childrenNodes[6 ^ a]));
                currNode = new_node(tx1, 8, ty1, 8, tzm, 7);
                break; }
            case 7: {
                nodeQueue.push(queueNode(txm, tym, tzm, tx1, ty1, tz1, node->childrenNodes[7 ^ a]));
                currNode = 8;
                break; }
            }
        } while (currNode < 8);
    }
}

void OCTREE::ray_parameter(OctreeNode* node, CRay ray) {
    // set tilt bits to 0
    a = 0;

    // all components of ray's direction vector must be positive
    // if input ray has some some negatives, their direction and
    // octree nodes will be tilted respectively

    // X direction component
    if (ray.GetDir().getX() < 0) {
        ray.SetLocX(-ray.GetLoc().getX());
        ray.SetDirX(-ray.GetDir().getX());
        a |= 4;
    }

    // Y direction component
    if (ray.GetDir().getY() < 0) {
        ray.SetLocY( -ray.GetLoc().getY());
        ray.SetDirY(-ray.GetDir().getY());
        a |= 2;
    }

    // Z direction component
    if (ray.GetDir().getZ() < 0) {
        ray.SetLocZ(-ray.GetLoc().getZ());
        ray.SetDirZ(-ray.GetDir().getZ());
        a |= 1;
    }
    
    float tx0 = (bbox.MinX() - ray.GetLoc().getX()) / ray.GetDir().getX();
    float tx1 = (bbox.MaxX() - ray.GetLoc().getX()) / ray.GetDir().getX();
    float ty0 = (bbox.MinY() - ray.GetLoc().getY()) / ray.GetDir().getY();
    float ty1 = (bbox.MaxY() - ray.GetLoc().getY()) / ray.GetDir().getY();
    float tz0 = (bbox.MinZ() - ray.GetLoc().getZ()) / ray.GetDir().getZ();
    float tz1 = (bbox.MaxZ() - ray.GetLoc().getZ()) / ray.GetDir().getZ();

    if (max(max(tx0, ty0), tz0) < min(min(tx1, ty1), tz1)) {
        proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, node);
        //proc_subtreeNonRec(tx0, ty0, tz0, tx1, ty1, tz1, node);
    }
}

const CObject3D* OCTREE::FindNearestI(CRay& ray, CHitPointInfo& info)
{
    numberOfQueries++;
    termNodesCnt = 0;
    endSearch = false;

    finalObj = nullptr;
    currentRay = &ray;
    currentInfo = &info;

    //info.maxt += CLimits::Threshold;
    //float saveMaxT = info->maxt;

    // calling the parametric octree algorithm
    ray_parameter(root, ray);

    //info.t = info.maxt;
    // recover maximum signed distance
    //info.maxt = saveMaxT;

    return (CObject3D*)(info.object = finalObj);

    //currentInfo->t = currentInfo->maxt;
    //currentInfo->maxt = saveMaxT;
}

void OCTREE::ProvideID(ostream& app)
{
    app << "OCTREE - fast recursive algorithm\n";
    app << (int)GetID() << "\n";
}

void OCTREE::printTriagle(CObject3D* object) {
    CTriangle* tmp = static_cast<CTriangle*>(object);
    std::cout << tmp->vertices[0].x << " " << tmp->vertices[0].y << " " << tmp->vertices[0].z << " " <<
        tmp->vertices[1].x << " " << tmp->vertices[1].y << " " << tmp->vertices[1].z << " " <<
        tmp->vertices[2].x << " " << tmp->vertices[2].y << " " << tmp->vertices[2].z << std::endl;
}

std::string printVertex(CVector3D& vertex)
{
    return to_string(vertex.x) + " " + to_string(vertex.y) + " " + to_string(vertex.z);
}

float roundOwn(float var)
{
    float value = (float)(int)(var * 100 + .5);
    return (float)value / 100;
};

void OCTREE::printObjects(const CObjectList& objlist) {
    int i = 0;
    std::cout << "--------------------------------------------------------" << std::endl;
    for (CObjectList::iterator it = ((CObjectList*)&objlist)->begin(); it != objlist.end(); it++, i++) {
        if (i >= 0)
            break;

        CObject3D* obj = *it;
        //objbox = obj->GetBox();
        std::cout << "Object: " << obj->UniqueID() << " is in " <<
            setfill(' ') << setw(7) << roundOwn(static_cast<CTriangle*>(obj)->GetCenter(0)) << " " <<
            setfill(' ') << setw(7) << roundOwn(static_cast<CTriangle*>(obj)->GetCenter(1)) << " " <<
            setfill(' ') << setw(7) << roundOwn(static_cast<CTriangle*>(obj)->GetCenter(2)) << std::endl;
    }
    std::cout << "--------------------------------------------------------" << std::endl;;
}

void OCTREE::printAllMeshes(std::vector<CObject3D*> objlist, string color) {

    for (auto& it : objlist) {
        CTriangle* triangleObj = static_cast<CTriangle*>(it);

        myfile << "glcolor " << color << "\n";
        myfile << "raw_triangle" << "\n";
        myfile << printVertex(triangleObj->vertices[0]) << " "
            << printVertex(triangleObj->vertices[1]) << " "
            << printVertex(triangleObj->vertices[2]) << "\n";
        myfile << "raw_end" << "\n";

    }
}

void OCTREE::printAllMeshes(const CObjectList& objlist) {

    for (CObjectList::iterator it = ((CObjectList*)&objlist)->begin(); it != objlist.end(); it++) {
        CTriangle* triangleObj = static_cast<CTriangle*>(*it);

        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        myfile << "glcolor " << to_string(r) << " " << to_string(r) << " " << to_string(r) << "\n";
        myfile << "raw_triangle" << "\n";
        myfile << printVertex(triangleObj->vertices[0]) << " "
            << printVertex(triangleObj->vertices[1]) << " "
            << printVertex(triangleObj->vertices[2]) << "\n";
        myfile << "raw_end" << "\n";
    }
}

void OCTREE::printRay(CRay& ray) {
    float t = bbox.Diagonal().MaxComponent() * 6;

    myfile << "glcolor 1 1 1" << "\n";
    myfile << "gllinewidth " << 7 << "\n";
    myfile << "raw_line" << "\n";

    myfile << printVertex(ray.GetLoc()) << " " << printVertex(ray.GetLoc() + ray.GetDir() * t) << "\n";

    myfile << "raw_end" << "\n";
}

void OCTREE::printTriangle(CObject3D& triangle, string color) {
    CTriangle* triangleObj = static_cast<CTriangle*>(&triangle);

    triangleObj->vertices[0];
    triangleObj->vertices[1];
    triangleObj->vertices[2];

    myfile << "glcolor " << color << "\n";
    myfile << "gllinewidth " << "1000" << "\n";
    myfile << "raw_triangle" << "\n";
    
    myfile << printVertex(triangleObj->vertices[0]) << " " 
           << printVertex(triangleObj->vertices[1]) << " " 
           << printVertex(triangleObj->vertices[2]) << "\n";

    myfile << "raw_end" << "\n";
}

void OCTREE::traverseOctree(OctreeNode* current) {
    if (current->terminal && GENERATE_MODEL) {
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float g = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float b = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        string color = to_string(r) + " " + to_string(g) + " " + to_string(b);

        //printVerticesFromBBox(current->nodeBBox, color, "1");
        printAllMeshes(*(current->objects), color);
    }

    for (size_t i = 0; i < 8; i++) {
        if (current->childrenNodes[i]) {
            traverseOctree(current->childrenNodes[i]);
        }
    }
}

void OCTREE::printVerticesFromBBox(const SBBox& boxPrint, string RGBcolor, string lineWidth) {
    CVector3D vertices[8];
    vertices[0] = boxPrint.Min();

    vertices[1] = CVector3D(boxPrint.MaxX(), boxPrint.MinY(), boxPrint.MinZ());
    vertices[2] = CVector3D(boxPrint.MinX(), boxPrint.MinY(), boxPrint.MaxZ());
    vertices[3] = CVector3D(boxPrint.MaxX(), boxPrint.MinY(), boxPrint.MaxZ());

    vertices[4] = CVector3D(boxPrint.MinX(), boxPrint.MaxY(), boxPrint.MinZ());
    vertices[5] = CVector3D(boxPrint.MaxX(), boxPrint.MaxY(), boxPrint.MinZ());
    vertices[6] = CVector3D(boxPrint.MinX(), boxPrint.MaxY(), boxPrint.MaxZ());

    vertices[7] = boxPrint.Max();

    myfile << "glcolor " + RGBcolor << "\n";
    myfile << "gllinewidth " << lineWidth << "\n";
    myfile << "raw_line" << "\n";

    myfile << printVertex(vertices[0]) << " " << printVertex(vertices[1]) << "\n";
    myfile << printVertex(vertices[0]) << " " << printVertex(vertices[2]) << "\n";
    myfile << printVertex(vertices[2]) << " " << printVertex(vertices[3]) << "\n";
    myfile << printVertex(vertices[1]) << " " << printVertex(vertices[3]) << "\n";

    myfile << printVertex(vertices[4]) << " " << printVertex(vertices[5]) << "\n";
    myfile << printVertex(vertices[4]) << " " << printVertex(vertices[6]) << "\n";
    myfile << printVertex(vertices[5]) << " " << printVertex(vertices[7]) << "\n";
    myfile << printVertex(vertices[6]) << " " << printVertex(vertices[7]) << "\n";

    myfile << printVertex(vertices[0]) << " " << printVertex(vertices[4]) << "\n";
    myfile << printVertex(vertices[1]) << " " << printVertex(vertices[5]) << "\n";
    myfile << printVertex(vertices[2]) << " " << printVertex(vertices[6]) << "\n";
    myfile << printVertex(vertices[3]) << " " << printVertex(vertices[7]) << "\n";

    myfile << "raw_end" << "\n";
}